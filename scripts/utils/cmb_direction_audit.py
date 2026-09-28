#!/usr/bin/env python3
"""
CMB directional-template audit for Step 2.7 (vectorized recompute).

Background
----------
The orbital-velocity term in calculate_3d_velocity_vectors was found to be
antiparallel to Earth's true heliocentric velocity (a 180-degree sign error):
at vernal equinox the coded vector pointed to ecliptic longitude 90 deg instead
of 270 deg. Because the predictor cos(dec of v_orb + v_bg) is invariant under
v_net -> -v_net, the buggy landscape is exactly the corrected landscape
mirrored through the origin: r_buggy(B) = r_correct(antipode B). This script
recomputes the per-combination best-fit directions, control directions and the
|V_bg| scan under the corrected sign, verifies the mirror relation, quantifies
the hemisphere degeneracy of the cos(dec) predictor, and writes an audit JSON.

Outputs:
    - results/outputs/step_2_7_cmb_direction_audit.json
"""

import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from datetime import date

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts" / "steps"))

OUTPUTS_DIR = ROOT / "results" / "outputs"
OUT_JSON = OUTPUTS_DIR / "step_2_7_cmb_direction_audit.json"

EPS = np.radians(23.44)
GRID_SPEED = 20.0
SPEED_SCAN = [5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 369.82]

CONTROL_DIRS = {
    'cmb_dipole':      (167.94, -6.94),
    'cmb_antipode':    (347.94, +6.94),
    'solar_apex':      (271.0, 30.0),
    'solar_antapex':   (91.0, -30.0),
    'ecliptic_np':     (270.0, 66.56),
    'ecliptic_sp':     (90.0, -66.56),
    'galactic_center': (266.4, -29.0),
    'perihelion_vel':  (191.87, -5.09),
    'aphelion_vel':    (11.87, +5.09),
}


def v_orb_ecl(doy, speed):
    """Earth heliocentric orbital velocity, ecliptic frame (CORRECTED sign).

    At vernal equinox (doy~80) Earth is at ecliptic longitude 180 deg and its
    prograde velocity points to lambda = 270 deg, i.e. (0, -1, 0) * speed.
    """
    th = (np.asarray(doy, float) - 80.0) / 365.25 * 2 * np.pi
    return np.stack([np.sin(th) * speed, -np.cos(th) * speed,
                     np.zeros_like(th)], axis=-1)


def ecl_to_eq(vec):
    x, y, z = vec[..., 0], vec[..., 1], vec[..., 2]
    return np.stack([x, y * np.cos(EPS) - z * np.sin(EPS),
                     y * np.sin(EPS) + z * np.cos(EPS)], axis=-1)


def net_dec(doy, speed, bra, bdec, bspeed, corrected=True):
    """Declination (deg) of v_orb + v_bg for each doy."""
    br, bd = np.radians(bra), np.radians(bdec)
    bg = np.array([bspeed * np.cos(bd) * np.cos(br),
                   bspeed * np.cos(bd) * np.sin(br),
                   bspeed * np.sin(bd)])
    vo = ecl_to_eq(v_orb_ecl(doy, np.asarray(speed, float)))
    if not corrected:
        vo = -vo
    v = vo + bg
    return np.degrees(np.arcsin(v[..., 2] / np.linalg.norm(v, axis=-1)))


def corr_at(doy, speed, ratio, bra, bdec, bspeed, corrected=True):
    p = np.cos(np.radians(net_dec(doy, speed, bra, bdec, bspeed, corrected)))
    if np.std(p) < 1e-12:
        return 0.0
    return float(stats.pearsonr(p, ratio)[0])


def best_direction(doy, speed, ratio, bspeed, ra_step=1.0, dec_step=1.0,
                   corrected=True):
    """Vectorized grid search of cos(dec_net) correlation."""
    ra_g = np.arange(0.0, 360.0, ra_step)
    dec_g = np.arange(-90.0, 90.0 + dec_step / 2, dec_step)
    RA, DEC = np.meshgrid(ra_g, dec_g)
    br, bd = np.radians(RA), np.radians(DEC)
    B = np.stack([bspeed * np.cos(bd) * np.cos(br),
                  bspeed * np.cos(bd) * np.sin(br),
                  bspeed * np.sin(bd)], axis=-1)          # (ndec, nra, 3)
    vo = ecl_to_eq(v_orb_ecl(doy, np.asarray(speed, float)))  # (nm, 3)
    if not corrected:
        vo = -vo
    V = vo[:, None, None, :] + B[None, :, :, :]             # (nm, ndec, nra, 3)
    pred = np.cos(np.arcsin(V[..., 2] / np.linalg.norm(V, axis=-1)))
    # Pearson r against observed ratio, vectorized
    p = pred.reshape(len(doy), -1)
    p = p - p.mean(axis=0)
    o = ratio - ratio.mean()
    denom = np.sqrt((p ** 2).sum(axis=0)) * np.sqrt((o ** 2).sum())
    r = np.where(denom > 0, (p * o[:, None]).sum(axis=0) / denom, 0.0)
    r = r.reshape(pred.shape[1:])
    i = np.nanargmax(r)
    di, ri = np.unravel_index(i, r.shape)
    return {'ra': float(ra_g[ri]), 'dec': float(dec_g[di]),
            'r': float(r[di, ri])}


def ang_sep(a, b, c, d):
    x = (np.sin(np.radians(b)) * np.sin(np.radians(d)) +
         np.cos(np.radians(b)) * np.cos(np.radians(d)) * np.cos(np.radians(a - c)))
    return float(np.degrees(np.arccos(np.clip(x, -1, 1))))


def load_monthly(filter_name='all_stations', mode='multi_gnss',
                 metric='pos_jitter', coherence='phase_alignment'):
    f = OUTPUTS_DIR / f"step_2_5_orbital_coupling_{filter_name}_{mode}.json"
    data = json.load(open(f))
    ma = data['results'][metric][coherence]['monthly_anisotropy']
    rows = []
    for k, v in sorted(ma.items()):
        rows.append((date(v['year'], v['month'], 15).timetuple().tm_yday,
                     v['ratio'], v.get('velocity', 29.78)))
    rows.sort()
    return (np.array([r[0] for r in rows], float),
            np.array([r[1] for r in rows], float),
            np.array([r[2] for r in rows], float))


def mirror_check(doy, speed, ratio, n=200):
    """Verify r_buggy(B) = r_correct(antipode B) at random directions."""
    rng = np.random.default_rng(0)
    maxdiff = 0.0
    for _ in range(n):
        ra = rng.uniform(0, 360)
        de = rng.uniform(-90, 90)
        bs = float(rng.choice(SPEED_SCAN))
        a = corr_at(doy, speed, ratio, ra, de, bs, corrected=False)
        b = corr_at(doy, speed, ratio, (ra + 180) % 360, -de, bs, corrected=True)
        maxdiff = max(maxdiff, abs(a - b))
    return maxdiff


def hemisphere_check(doy, speed, ratio, bspeed=20.0, corrected=True):
    """cos(dec_net) is even in net declination: compare r(B) for the bg
    directions mirrored about the equatorial plane at the corrected best-fit
    region. Report the maximum |r| difference over a grid — a measure of how
    weakly the hemisphere of the fixed vector is constrained."""
    ra_g = np.arange(0.0, 360.0, 10.0)
    dec_g = np.arange(-80.0, 80.0 + 1, 10.0)
    diffs, same_sign = [], 0
    for de in dec_g:
        if de == 0:
            continue
        for ra in ra_g:
            rn = corr_at(doy, speed, ratio, ra, abs(de), bspeed, corrected)
            rs = corr_at(doy, speed, ratio, ra, -abs(de), bspeed, corrected)
            diffs.append(abs(rn - rs))
            if rn * rs > 0:
                same_sign += 1
    return {'max_abs_r_diff_ns_mirror': float(np.max(diffs)),
            'median_abs_r_diff_ns_mirror': float(np.median(diffs)),
            'frac_same_sign': float(same_sign / len(diffs))}


def main():
    audit = {'step': '2.7-audit',
             'name': 'cmb_direction_template_audit',
             'bug': ('Orbital velocity vector in calculate_3d_velocity_vectors '
                     'was antiparallel to the true Earth velocity (180-deg '
                     'sign error: coded v_ecl=(0,+1,0) at vernal equinox vs '
                     'correct (0,-1,0)). Fixed in step_2_7; this audit '
                     'recomputes all directional claims.'),
             'mirror_identity': ('Predictor cos(dec(v_net)) is invariant under '
                                 'v_net -> -v_net, so r_buggy(B) = '
                                 'r_correct(antipode(B)) exactly.')}

    # primary combination (the one used by step_2_7_controls_speed_scan.json)
    doy, ratio, spd = load_monthly()
    audit['data'] = {'n_months': len(doy), 'span_years': 3,
                     'combo': 'all_stations/multi_gnss/pos_jitter/phase_alignment'}

    audit['mirror_check_max_abs_diff'] = mirror_check(doy, spd, ratio)

    # corrected landscape: fine grid at 20 km/s
    best = best_direction(doy, spd, ratio, GRID_SPEED)
    best['ecliptic_lat_deg'] = None
    # ecliptic latitude of corrected best
    ra_, de_ = np.radians(best['ra']), np.radians(best['dec'])
    sinb = (np.sin(de_) * np.cos(EPS)
            - np.cos(de_) * np.sin(EPS) * np.sin(ra_))
    best['ecliptic_lat_deg'] = float(np.degrees(np.arcsin(np.clip(sinb, -1, 1))))
    best['sep_cmb_deg'] = ang_sep(best['ra'], best['dec'], 167.94, -6.94)
    best['sep_aphelion_tangent_deg'] = ang_sep(best['ra'], best['dec'],
                                             11.87, 5.09)
    best['sep_perihelion_tangent_deg'] = ang_sep(best['ra'], best['dec'],
                                               191.87, -5.09)
    audit['corrected_best_fit_20kms'] = best

    # controls under corrected sign
    ctrl = {}
    for name, (ra, de) in CONTROL_DIRS.items():
        ctrl[name] = {'ra': ra, 'dec': de,
                      'r_corrected': corr_at(doy, spd, ratio, ra, de,
                                             GRID_SPEED, corrected=True),
                      'r_buggy': corr_at(doy, spd, ratio, ra, de,
                                         GRID_SPEED, corrected=False)}
    audit['control_directions_20kms'] = ctrl

    # speed scan under corrected sign (2-deg grid for speed)
    scan = []
    for bs in SPEED_SCAN:
        b = best_direction(doy, spd, ratio, bs, ra_step=1.0, dec_step=1.0)
        b['speed_kms'] = bs
        b['cmb_sep_deg'] = ang_sep(b['ra'], b['dec'], 167.94, -6.94)
        b['aphelion_tangent_sep_deg'] = ang_sep(b['ra'], b['dec'], 11.87, 5.09)
        scan.append(b)
    audit['corrected_speed_scan'] = scan

    # hemisphere degeneracy of the predictor
    audit['hemisphere_degeneracy'] = hemisphere_check(doy, spd, ratio)

    # corrected landscape on every available combination (coarse 5-deg grid)
    combos = {}
    for filt in ['all_stations', 'optimal_100', 'dynamic_50']:
        for mode in ['baseline', 'multi_gnss', 'ionofree']:
            for metric in ['clock_bias', 'pos_jitter', 'clock_drift']:
                for coh in ['msc', 'phase_alignment']:
                    f = OUTPUTS_DIR / f"step_2_5_orbital_coupling_{filt}_{mode}.json"
                    if not f.exists():
                        continue
                    data = json.load(open(f))
                    try:
                        ma = data['results'][metric][coh]['monthly_anisotropy']
                    except KeyError:
                        continue
                    rows = []
                    for k, v in sorted(ma.items()):
                        if 'ratio' not in v:
                            continue
                        rows.append((date(v['year'], v['month'], 15)
                                     .timetuple().tm_yday, v['ratio'],
                                     v.get('velocity', 29.78)))
                    if len(rows) < 12:
                        continue
                    rows.sort()
                    d_ = np.array([r[0] for r in rows], float)
                    r_ = np.array([r[1] for r in rows], float)
                    s_ = np.array([r[2] for r in rows], float)
                    bc = best_direction(d_, s_, r_, GRID_SPEED,
                                        ra_step=5.0, dec_step=5.0)
                    bb = best_direction(d_, s_, r_, GRID_SPEED,
                                        ra_step=5.0, dec_step=5.0,
                                        corrected=False)
                    key = f"{filt}/{mode}/{metric}/{coh}"
                    combos[key] = {
                        'corrected': {**bc,
                                      'cmb_sep_deg': ang_sep(
                                          bc['ra'], bc['dec'], 167.94, -6.94)},
                        'buggy': {**bb,
                                  'cmb_sep_deg': ang_sep(
                                      bb['ra'], bb['dec'], 167.94, -6.94)},
                        'n_months': len(d_),
                    }
    audit['all_combinations'] = combos
    n_cmb_45 = sum(1 for c in combos.values()
                   if c['corrected']['cmb_sep_deg'] < 45)
    n_aph_10 = sum(1 for c in combos.values()
                   if ang_sep(c['corrected']['ra'], c['corrected']['dec'],
                              11.87, 5.09) < 15)
    audit['combination_summary'] = {
        'n_combinations': len(combos),
        'n_corrected_within_45deg_of_cmb': n_cmb_45,
        'n_corrected_within_15deg_of_aphelion_tangent': n_aph_10,
    }

    OUT_JSON.write_text(json.dumps(audit, indent=2))
    print(json.dumps({'best': best,
                      'mirror_diff': audit['mirror_check_max_abs_diff'],
                      'hemi': audit['hemisphere_degeneracy'],
                      'combos': audit['combination_summary']}, indent=2))


if __name__ == '__main__':
    main()
