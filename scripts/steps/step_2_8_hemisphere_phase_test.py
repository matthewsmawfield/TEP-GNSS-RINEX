#!/usr/bin/env python3
"""step_2_8_hemisphere_phase_test: Hemisphere-stratified seasonal-phase test and harmonic decomposition.

Issue 3-2 diagnostic: distinguishes a heliocentric (orbital-velocity) driver of
the monthly anisotropy modulation from a local-seasonal (insolation/ionosphere)
driver. A heliocentric driver predicts the same annual phase in both
hemispheres; a local-seasonal driver predicts anti-phase (~pi) between them.

Additionally decomposes the monthly series into annual (A1) and semiannual (A2)
harmonics to discriminate between orbital velocity coupling (pure annual harmonic,
perihelion phase-locked) and equinoctial systematics (Russell-McPherron geomagnetic
coupling and GPS eclipse seasons).

Inputs (existing Step 2.5 outputs, no refitting):
  results/outputs/step_2_5_orbital_coupling_dynamic_50_north_multi_gnss.json
  results/outputs/step_2_5_orbital_coupling_dynamic_50_south_multi_gnss.json
  results/outputs/step_2_5_orbital_coupling_all_stations_multi_gnss.json
  results/outputs/step_2_5_orbital_coupling_dynamic_50_multi_gnss.json
  results/outputs/step_2_5_orbital_coupling_dynamic_50_baseline.json
  results/outputs/step_2_5_orbital_coupling_dynamic_50_ionofree.json

Output: results/outputs/step_2_8_hemisphere_phase_test.json
"""
import hashlib
import json
import math
import os
import sys
from datetime import date
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import eval_legendre

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0,ROOT)
from scripts.steps.step_2_0_raw_spp_analysis import _compute_coherence_phase,FS_HZ,F1_HZ,F2_HZ
OUT = os.path.join(ROOT, "results", "outputs")
NORTH = os.path.join(OUT, "step_2_5_orbital_coupling_dynamic_50_north_multi_gnss.json")
SOUTH = os.path.join(OUT, "step_2_5_orbital_coupling_dynamic_50_south_multi_gnss.json")

METRICS = ["clock_bias", "pos_jitter"]
SUBS = ["msc", "phase_alignment"]


def load_series(path, metric, sub):
    d = json.load(open(path))
    m = d["results"][metric][sub]["monthly_anisotropy"]
    out = {}
    for k, v in m.items():
        if not isinstance(v, dict):
            continue
        r = v.get("ratio") or v.get("anisotropy_ratio")
        if r is None:
            ew, ns = v.get("ew_lambda"), v.get("ns_lambda")
            if ew and ns:
                r = ew / ns
        if r is not None and math.isfinite(r):
            out[k] = float(r)
    return out


def fit_annual(months, y):
    """OLS fit y ~ 1 + sin(2 pi t / 12) + cos(2 pi t / 12), t in months."""
    t = np.array([(int(k[:4]) - 2022) * 12 + int(k[5:7]) - 1 for k in months], float)
    y = np.asarray(y, float)
    X = np.column_stack([np.ones_like(t), np.sin(2 * np.pi * t / 12.0), np.cos(2 * np.pi * t / 12.0)])
    beta, res, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = max(len(y) - 3, 1)
    sigma2 = float(resid @ resid / dof)
    cov = sigma2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    amp = float(math.hypot(beta[1], beta[2]))
    phase = float(math.atan2(beta[1], beta[2]))  # atan2(sin_coef, cos_coef)
    phase_deg = (math.degrees(phase) + 360.0) % 360.0
    amp_se = float(math.sqrt(max(se[1] ** 2 * math.cos(phase) ** 2 + se[2] ** 2 * math.sin(phase) ** 2, 0.0))) if amp else float("nan")
    rss1 = float(resid @ resid)
    rss0 = float(((y - y.mean()) ** 2).sum())
    df1, df2 = 2, dof
    F = ((rss0 - rss1) / df1) / (rss1 / df2) if rss1 > 0 else float("inf")
    try:
        from scipy.stats import f as _f
        p = float(_f.sf(F, df1, df2))
    except Exception:
        p = float("nan")
    return {
        "amplitude": amp,
        "amplitude_se": amp_se,
        "phase_rad": phase,
        "phase_deg": phase_deg,
        "F": float(F),
        "p_value": p,
        "const": float(beta[0]),
    }


def fit_harmonics(months, y):
    """Decompose y(t) into constant, annual harmonic (A1), and semiannual harmonic (A2)."""
    t = np.array([(int(k[:4]) - 2022) * 12 + int(k[5:7]) - 1 for k in months], float)
    y = np.asarray(y, float)
    w = 2.0 * np.pi / 12.0

    X_ann = np.column_stack([np.ones_like(t), np.cos(w * t), np.sin(w * t)])
    X_full = np.column_stack([np.ones_like(t), np.cos(w * t), np.sin(w * t), np.cos(2 * w * t), np.sin(2 * w * t)])

    b_ann, *_ = np.linalg.lstsq(X_ann, y, rcond=None)
    b_full, *_ = np.linalg.lstsq(X_full, y, rcond=None)

    rss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2_ann = float(1.0 - np.sum((y - X_ann @ b_ann) ** 2) / rss_tot) if rss_tot > 0 else 0.0
    r2_full = float(1.0 - np.sum((y - X_full @ b_full) ** 2) / rss_tot) if rss_tot > 0 else 0.0

    c0 = float(b_full[0])
    A1 = float(math.hypot(b_full[1], b_full[2]))
    phi1 = float(math.atan2(b_full[2], b_full[1]))
    A2 = float(math.hypot(b_full[3], b_full[4]))
    phi2 = float(math.atan2(b_full[4], b_full[3]))

    # Annual peak and trough timing in 0-indexed month within year [0=Jan, 1=Feb, ...]
    t_peak = float((phi1 / w) % 12)
    t_trough = float(((phi1 + math.pi) / w) % 12)
    doy_trough = float((t_trough * 30.4375 + 1) % 365.25)
    doy_peak = float((t_peak * 30.4375 + 1) % 365.25)

    # Offset from Perihelion (DOY 4.0, Jan 4)
    diff_peri = (4.0 - doy_trough + 365.25) % 365.25
    if diff_peri > 182.6:
        diff_peri -= 365.25

    # Offset from Winter Solstice (DOY 355.0, Dec 21)
    diff_sol = (355.0 - doy_trough + 365.25) % 365.25
    if diff_sol > 182.6:
        diff_sol -= 365.25

    # Semiannual peaks
    p1_m = float((phi2 / (2 * w)) % 6)
    p2_m = float((p1_m + 6) % 12)
    doy_semi1 = float((p1_m * 30.4375 + 1) % 365.25)
    doy_semi2 = float((p2_m * 30.4375 + 1) % 365.25)

    # Offset from Vernal Equinox (DOY 80.0, March 21) and Autumnal Equinox (DOY 265.0, Sept 22)
    diff_vernal = abs(doy_semi1 - 80.0) if doy_semi1 < 180 else abs(doy_semi2 - 80.0)
    diff_autumn = abs(doy_semi2 - 265.0) if doy_semi2 >= 180 else abs(doy_semi1 - 265.0)

    ratio_a1_a2 = float(A1 / max(A2, 1e-6))

    return {
        "c0": c0,
        "A1_annual": A1,
        "A2_semiannual": A2,
        "A1_over_A2": ratio_a1_a2,
        "R2_annual": r2_ann,
        "R2_annual_plus_semiannual": r2_full,
        "annual_trough_month": t_trough,
        "annual_trough_doy": doy_trough,
        "annual_peak_month": t_peak,
        "annual_peak_doy": doy_peak,
        "offset_perihelion_days": float(diff_peri),
        "offset_perihelion_deg": float((diff_peri / 365.25) * 360.0),
        "offset_solstice_days": float(diff_sol),
        "offset_solstice_deg": float((diff_sol / 365.25) * 360.0),
        "semiannual_peak1_doy": doy_semi1,
        "semiannual_peak2_doy": doy_semi2,
        "offset_vernal_equinox_days": float(diff_vernal),
        "offset_autumnal_equinox_days": float(diff_autumn),
    }


def wrap_pi(dphi):
    return (dphi + math.pi) % (2 * math.pi) - math.pi


def locked_orbital_replay(path, metric, sub):
    entries = json.load(open(path))["results"][metric][sub]["monthly_anisotropy"]
    rows = {}
    for record in entries.values():
        if not isinstance(record, dict):
            continue
        ratio, velocity = record.get("ratio"), record.get("velocity")
        if ratio and velocity and ratio > 0 and math.isfinite(ratio) and math.isfinite(velocity):
            rows[(int(record["year"]), int(record["month"]))] = (
                float(ratio), float(velocity),
                int(record.get("n_ew_pairs",0)), int(record.get("n_ns_pairs",0)))
    if any((year, month) not in rows for year in (2022, 2023, 2024) for month in range(1, 13)):
        return {"status": "incomplete", "months_available": len(rows)}

    def year_data(year, shift=0):
        y = np.array([rows[(year, month)][0] for month in range(1, 13)])
        v = np.array([rows[(year, (month+shift-1) % 12+1)][1]
                      for month in range(1, 13)])
        return y, v

    def evaluate(shift):
        y_train, v_train = year_data(2022, shift)
        center = float(v_train.mean())
        fit = np.linalg.lstsq(np.column_stack((np.ones(12), v_train-center)),
                              np.log(y_train), rcond=None)[0]
        annual = np.array([rows[(2022, month)][0] for month in range(1, 13)])
        evaluations = {}
        for year in (2023, 2024):
            observed, v = year_data(year, shift)
            predicted = np.exp(fit[0]+fit[1]*(v-center))
            flat = np.full(12, np.exp(fit[0]))
            same_year_slope = np.linalg.lstsq(
                np.column_stack((np.ones(12),v-v.mean())),
                np.log(observed),rcond=None)[0][1]
            coverage = np.array([rows[(year,month)][2]/max(rows[(year,month)][3],1)
                                 for month in range(1,13)])
            evaluations[str(year)] = {
                "observed_mean_ratio":float(np.mean(observed)),
                "mean_ew_to_ns_pair_count_ratio":float(np.mean(coverage)),
                "correlation_to_pair_coverage":float(np.corrcoef(observed,coverage)[0,1]) if np.std(coverage)>1e-12 else None,
                "descriptive_same_year_slope_per_km_s":float(same_year_slope),
                "rmse_velocity": float(np.sqrt(np.mean((observed-predicted)**2))),
                "rmse_constant": float(np.sqrt(np.mean((observed-flat)**2))),
                "rmse_calendar_twin": float(np.sqrt(np.mean((observed-annual)**2))),
                "correlation": float(np.corrcoef(observed,predicted)[0,1]) if (
                    np.std(predicted)>1e-12 and np.std(observed)>1e-12) else None,
                "observed_ew_gt_ns_months": int(np.sum(observed>1)),
                "predicted_ew_gt_ns_months": int(np.sum(predicted>1)),
                "n_months": 12,
            }
        return fit, evaluations

    fit, evaluations = evaluate(0)
    train_coverage = np.mean([rows[(2022,month)][2]/max(rows[(2022,month)][3],1)
                              for month in range(1,13)])
    shifts = [evaluate(shift)[1]["2024"]["rmse_velocity"] for shift in range(12)]
    return {"status":"retrospective internal replay, not blind or independent replication",
            "training_pair_count_ratio_ew_to_ns":float(train_coverage),
            "training_year":2022, "validation_year":2023, "test_year":2024,
            "model":"log(EW_lambda/NS_lambda) = intercept + slope * (orbital_speed - 2022 mean speed)",
            "fitted_intercept":float(fit[0]), "fitted_slope_per_km_s":float(fit[1]),
            "required_negative_slope":bool(fit[1]<0),
            "years":evaluations,
            "2024_velocity_phase_shift_rmse":shifts,
            "unshifted_phase_rank":int(1+sum(score<shifts[0] for score in shifts[1:])),
            "caveat":"Calendar-twin, insolation and orbital speed share the annual phase; a stable speed regression alone cannot identify a scalar source."}


def pair_level_kernel_replay():
    directory = Path(ROOT)/"data"/"processed"
    field_path = Path(ROOT).parent/"TEP"/"results"/"step_04_gamma_derivation.json"
    coords = json.load(open(directory/"station_coordinates.json"))
    field = json.load(open(field_path))["kinetic_forward_test"]
    stations = sorted(coords)
    xyz = np.asarray([coords[s] for s in stations],float)
    radii = np.linalg.norm(xyz,axis=1)
    valid = np.all(np.isfinite(xyz),axis=1)&(radii>6e6)&(radii<7e6)
    excluded_coordinates = int(np.sum(~valid))
    stations = [station for station,keep in zip(stations,valid) if keep]
    xyz = xyz[valid]/radii[valid,None]
    if not np.all(np.isfinite(xyz)) or np.max(np.abs(xyz))>1.01:
        raise ValueError("Station ECEF normalization failed")
    if len(stations)<40:
        return {"status":"fewer than 40 valid terrestrial ECEF stations",
                "valid_stations":len(stations)}
    selected = [0]
    spread = np.full(len(stations),np.inf)
    for _ in range(39):
        spread = np.minimum(spread,1-np.einsum("ij,j->i",xyz,xyz[selected[-1]]))
        spread[selected] = -1
        selected.append(int(np.argmax(spread)))
    stations = [stations[i] for i in selected]
    xyz = xyz[selected]
    angles = np.arccos(np.clip(np.einsum("ik,jk->ij",xyz,xyz),-1,1))
    source_cases = field["radial_source_scan"]
    if not all("angular_powers" in source_cases[name]
               for name in ("volume_white","core_only","mantle_only","surface_concentrated")):
        raise ValueError("Rebuild TEP step_04 before the source-family replay")
    powers = {name:np.asarray(case["angular_powers"],float)
              for name,case in source_cases.items()}
    ell = np.arange(1,len(powers["volume_white"]))
    weights = {name:(2*ell+1)*value[1:] for name,value in powers.items()}
    def kernel(theta,name="volume_white"):
        angular = [eval_legendre(int(l),math.cos(theta)) for l in ell]
        return float(np.dot(weights[name],angular)/weights[name].sum())
    bins = np.geomspace(50,13000,21)
    rows = {}
    for year in (2022,2023,2024):
        distances, values, months = [], [], []
        days_loaded = 0
        for month in range(1,13):
            doy = date(year,month,15).timetuple().tm_yday
            records = {}
            for station in stations:
                path = directory/f"{station}_{year}{doy:03d}.npz"
                if not path.exists():
                    continue
                with np.load(path) as product:
                    t = product.get("timestamps")
                    v = product.get("clock_bias_ns")
                    if t is None or v is None or len(t)!=len(v):
                        continue
                    valid = (t>=0)&(t<86400)&np.isfinite(v)
                    if np.sum(valid)>=200:
                        records[station]=(t[valid],v[valid])
            if len(records)<2:
                continue
            days_loaded += 1
            for i,sta1 in enumerate(stations):
                if sta1 not in records:
                    continue
                for j in range(i+1,len(stations)):
                    sta2 = stations[j]
                    distance = 6371*angles[i,j]
                    if sta2 not in records or not 50<=distance<=13000:
                        continue
                    t1,v1 = records[sta1]
                    t2,v2 = records[sta2]
                    _,a,b = np.intersect1d(t1,t2,return_indices=True)
                    if len(a)<200:
                        continue
                    _,_,phase = _compute_coherence_phase(v1[a],v2[b],FS_HZ,F1_HZ,F2_HZ)
                    if math.isfinite(phase):
                        distances.append(distance)
                        values.append(phase)
                        months.append(month)
        distance = np.asarray(distances)
        value = np.asarray(values)
        months = np.asarray(months)
        indices = np.searchsorted(bins,distance,side="right")-1
        summary = {}
        for index in range(len(bins)-1):
            mask = indices==index
            if np.sum(mask)>=6:
                summary[index]=(float(distance[mask].mean()),float(value[mask].mean()),int(np.sum(mask)))
        rows[year]={"bins":summary,"n_pair_days":len(value),"n_sample_days":days_loaded,
                    "distance":distance,"value":value,"month":months,"bin_index":indices}
    train = rows[2022]["bins"]
    if len(train)<5:
        return {"status":"insufficient long-baseline bins", "counts":{
            str(year):{key:rows[year][key] for key in ("n_pair_days","n_sample_days")}
            for year in rows}}
    r = np.array([train[i][0] for i in sorted(train)])
    y = np.array([train[i][1] for i in sorted(train)])
    angular_design = np.array([[eval_legendre(int(l),math.cos(d/6371))
                                for l in ell] for d in r])
    angular_rank = int(np.linalg.matrix_rank(angular_design))
    source_fits = {}
    source_training = {}
    for name in ("volume_white","core_only","mantle_only","surface_concentrated"):
        template = np.array([kernel(d/6371,name) for d in r])
        estimate = np.linalg.lstsq(np.column_stack((np.ones(len(r)),template)),y,rcond=None)[0]
        source_fits[name] = estimate
        source_training[name] = {"offset":float(estimate[0]),
            "amplitude":float(estimate[1]),
            "training_rmse":float(np.sqrt(np.mean((y-estimate[0]-estimate[1]*template)**2)))}
    fit = source_fits["volume_white"]
    linear = np.linalg.lstsq(np.column_stack((np.ones(len(r)),r)),y,rcond=None)[0]
    fixed_exp = np.linalg.lstsq(np.column_stack((np.ones(len(r)),
                                  np.exp(-r/4200))),y,rcond=None)[0]
    exponential = lambda distance,offset,amplitude,length: offset+amplitude*np.exp(-distance/length)
    free_exp,_ = curve_fit(exponential,r,y,p0=(y.mean(),.1,4200),
                           bounds=([-2,-2,100],[2,2,50000]),maxfev=10000)
    scores = {}
    for year in (2023,2024):
        shared = sorted(set(train)&set(rows[year]["bins"]))
        distance = np.array([rows[year]["bins"][i][0] for i in shared])
        observed = np.array([rows[year]["bins"][i][1] for i in shared])
        predicted = fit[0]+fit[1]*np.array([kernel(d/6371) for d in distance])
        source_scores = {}
        for name, estimate in source_fits.items():
            candidate = estimate[0]+estimate[1]*np.array([kernel(d/6371,name) for d in distance])
            source_scores[name] = float(np.sqrt(np.mean((observed-candidate)**2)))
        month_delete_improvement = []
        for excluded_month in range(1,13):
            keep = rows[year]["month"] != excluded_month
            selected_bins = rows[year]["bin_index"][keep]
            selected_r = rows[year]["distance"][keep]
            selected_y = rows[year]["value"][keep]
            means = [(float(selected_r[selected_bins==i].mean()),
                      float(selected_y[selected_bins==i].mean()))
                     for i in shared if np.sum(selected_bins==i)>=6]
            if len(means)<3:
                continue
            rr, yy = np.asarray(means).T
            kernel_pred = fit[0]+fit[1]*np.array([kernel(d/6371) for d in rr])
            linear_pred = linear[0]+linear[1]*rr
            month_delete_improvement.append(float(
                np.sqrt(np.mean((yy-linear_pred)**2))-
                np.sqrt(np.mean((yy-kernel_pred)**2))))
        scores[str(year)]={"n_bins":len(shared),
            "month_delete_improvement_linear_minus_kinetic":month_delete_improvement,
            "n_pair_days":rows[year]["n_pair_days"],
            "n_sample_days":rows[year]["n_sample_days"],
            "rmse_kinetic_shape":float(np.sqrt(np.mean((observed-predicted)**2))),
            "source_family_rmse":source_scores,
            "rmse_constant":float(np.sqrt(np.mean((observed-y.mean())**2))),
            "rmse_linear_distance":float(np.sqrt(np.mean((observed-linear[0]-linear[1]*distance)**2))),
            "rmse_fixed_4200km_exponential":float(np.sqrt(np.mean((observed-fixed_exp[0]-fixed_exp[1]*np.exp(-distance/4200))**2))),
            "rmse_train_fitted_exponential":float(np.sqrt(np.mean((observed-exponential(distance,*free_exp))**2))),
            "fitted_on_2022_amplitude":float(fit[1])}
    return {"status":"retrospective sparse-day same-network estimator-shape pilot, not a source-spectrum inversion or blind test",
        "station_selection":"40 deterministic farthest-point ECEF stations from sorted IDs with finite terrestrial radii, independent of clock data",
        "excluded_invalid_coordinates":excluded_coordinates,
        "station_ids":stations,
        "days":"15th of every month in 2022–2024",
        "processing":"baseline SPP clock_bias_ns, unchanged step_2_0 phase-alignment estimator",
        "field_kernel_source":"../TEP/results/step_04_gamma_derivation.json",
        "field_kernel_sha256":hashlib.sha256(field_path.read_bytes()).hexdigest(),
        "training_2022":{"n_bins":len(train),"n_pair_days":rows[2022]["n_pair_days"],
                         "n_sample_days":rows[2022]["n_sample_days"],
                         "bin_mean_distance_km":[train[i][0] for i in sorted(train)],
                         "offset":float(fit[0]),"amplitude":float(fit[1]),
                         "generic_exponential_fitted_length_km":float(free_exp[2])},
        "source_family_training_2022":source_training,
        "angular_source_identifiability":{
            "n_independent_training_bins":len(r),
            "n_angular_modes":len(ell),
            "matrix_rank_upper_bound":angular_rank,
            "unconstrained_null_directions_at_least":len(ell)-angular_rank,
            "scope":"Unregularized Legendre basis before nonnegative-power constraints and offset/transfer degeneracies; source spectrum not uniquely invertible from this pilot."},
        "validation":scores,
        "exploratory_controls":"The four radial source families were declared in TEP step_04; each offset and amplitude is trained only on 2022. Generic exponential length is trained only on 2022; fixed 4200 km is an existing GNSS calibration, not independent physical input.",
        "caveat":"One day per month, shared stations, fitted offset/amplitude and no processing-chain injection; cannot identify the physical source or claim an absolute covariance amplitude."}


def main():
    # 1. Hemisphere phase test
    north = json.load(open(NORTH))
    south = json.load(open(SOUTH))
    channels = {}
    for metric in METRICS:
        for sub in SUBS:
            N = load_series(NORTH, metric, sub)
            S = load_series(SOUTH, metric, sub)
            common = sorted(set(N) & set(S))
            if len(common) < 14:
                continue
            nv = [N[k] for k in common]
            sv = [S[k] for k in common]
            fN = fit_annual(common, nv)
            fS = fit_annual(common, sv)
            dphi = wrap_pi(fN["phase_rad"] - fS["phase_rad"])
            corr = float(np.corrcoef(nv, sv)[0, 1])
            channels[f"{metric}/{sub}"] = {
                "n_months": len(common),
                "north": fN,
                "south": fS,
                "phase_difference_rad": dphi,
                "phase_difference_deg": math.degrees(dphi),
                "monthly_correlation": corr,
                "in_phase": abs(dphi) < math.pi / 2,
            }
    n_in = sum(1 for c in channels.values() if c["in_phase"])

    # 2. Harmonic decomposition across representative configurations
    decomp_targets = [
        ("ALL_STATIONS_multi_gnss_pos_jitter_phase_alignment", os.path.join(OUT, "step_2_5_orbital_coupling_all_stations_multi_gnss.json"), "pos_jitter", "phase_alignment"),
        ("DYNAMIC_50_multi_gnss_pos_jitter_phase_alignment", os.path.join(OUT, "step_2_5_orbital_coupling_dynamic_50_multi_gnss.json"), "pos_jitter", "phase_alignment"),
        ("DYNAMIC_50_multi_gnss_clock_bias_msc", os.path.join(OUT, "step_2_5_orbital_coupling_dynamic_50_multi_gnss.json"), "clock_bias", "msc"),
        ("DYNAMIC_50_baseline_clock_bias_msc", os.path.join(OUT, "step_2_5_orbital_coupling_dynamic_50_baseline.json"), "clock_bias", "msc"),
        ("DYNAMIC_50_ionofree_clock_bias_msc", os.path.join(OUT, "step_2_5_orbital_coupling_dynamic_50_ionofree.json"), "clock_bias", "msc"),
    ]
    harmonic_results = {}
    for name, path, met, sub in decomp_targets:
        if os.path.exists(path):
            s = load_series(path, met, sub)
            months = sorted(s.keys())
            y = [s[k] for k in months]
            harmonic_results[name] = fit_harmonics(months, y)

    replay_targets = [
        ("primary_all_multi_clock_phase", "all_stations_multi_gnss", "clock_bias", "phase_alignment"),
        ("all_multi_position_phase", "all_stations_multi_gnss", "pos_jitter", "phase_alignment"),
        ("dynamic_multi_clock_msc", "dynamic_50_multi_gnss", "clock_bias", "msc"),
        ("dynamic_baseline_clock_msc", "dynamic_50_baseline", "clock_bias", "msc"),
        ("dynamic_ionofree_clock_msc", "dynamic_50_ionofree", "clock_bias", "msc"),
    ]
    locked_replay = {
        name:locked_orbital_replay(
            os.path.join(OUT, f"step_2_5_orbital_coupling_{product}.json"),metric,sub)
        for name, product, metric, sub in replay_targets
    }

    result = {
        "step": "2_8_hemisphere_phase_test",
        "purpose": "Distinguish heliocentric (in-phase) vs local-seasonal (anti-phase) driver, and decompose annual vs semiannual power",
        "inputs": {
            "north": os.path.basename(NORTH),
            "south": os.path.basename(SOUTH),
            "filter": "DYNAMIC_50",
            "mode": "multi_gnss",
        },
        "channels": channels,
        "harmonic_decomposition": harmonic_results,
        "locked_orbital_replay":locked_replay,
        "summary": {
            "n_channels": len(channels),
            "n_in_phase": n_in,
            "locked_primary_2024_outcome": "fails fixed-template prediction" if (
                locked_replay["primary_all_multi_clock_phase"]["years"]["2024"]["rmse_velocity"] >=
                min(locked_replay["primary_all_multi_clock_phase"]["years"]["2024"]["rmse_constant"],
                    locked_replay["primary_all_multi_clock_phase"]["years"]["2024"]["rmse_calendar_twin"])
                or (locked_replay["primary_all_multi_clock_phase"]["years"]["2024"]["correlation"] or 0) <= 0
            ) else "passes these internal comparators",
            "verdict": (
                f"{n_in}/{len(channels)} channels have full-sample annual phase differences below 90 degrees; "
                "the Multi-GNSS position-phase series has a strong annual component, while clock MSC retains semiannual power. "
                "These are descriptive full-sample results, not a source identification: shared seasonal and processing effects "
                "can remain in phase across hemispheres. The separately reported fixed 2022 orbital-speed regression must pass "
                "the 2023 and 2024 constant, calendar-twin and phase-shift controls before predictive orbital support is claimed."
            ),
        },
    }
    if "--kernel-pilot" in sys.argv:
        result["pair_level_kernel_replay"] = pair_level_kernel_replay()
    out_path = os.path.join(OUT, "step_2_8_hemisphere_phase_test.json")
    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=2)

    print("=== Hemisphere Phase Test ===")
    for name, c in channels.items():
        print(
            f"{name}: N amp={c['north']['amplitude']:.3f} ph={c['north']['phase_deg']:.0f}deg | "
            f"S amp={c['south']['amplitude']:.3f} ph={c['south']['phase_deg']:.0f}deg | "
            f"dphi={c['phase_difference_deg']:+.0f}deg corr={c['monthly_correlation']:+.2f} "
            f"{'IN-PHASE' if c['in_phase'] else 'ANTI-PHASE/AMBIGUOUS'}"
        )
    print("\n=== Harmonic Decomposition (Annual vs Semiannual) ===")
    for name, h in harmonic_results.items():
        print(
            f"{name}:\n"
            f"  A1={h['A1_annual']:.3f}, A2={h['A2_semiannual']:.3f}, A1/A2={h['A1_over_A2']:.2f}, R2_ann={h['R2_annual']:.2f}, R2_full={h['R2_annual_plus_semiannual']:.2f}\n"
            f"  Trough DOY={h['annual_trough_doy']:.1f} (offset Perihelion: {h['offset_perihelion_days']:+.1f}d, offset Solstice: {h['offset_solstice_days']:+.1f}d)\n"
            f"  Semiannual peaks: DOY {h['semiannual_peak1_doy']:.1f}, {h['semiannual_peak2_doy']:.1f} (Vernal d={h['offset_vernal_equinox_days']:.1f}d, Autumn d={h['offset_autumnal_equinox_days']:.1f}d)"
        )
    print("\n" + result["summary"]["verdict"])
    print("wrote", out_path)


if __name__ == "__main__":
    main()
