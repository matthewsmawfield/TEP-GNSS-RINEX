#!/usr/bin/env python3
"""
TEP-GNSS-RINEX Analysis - STEP 2.6b: Season-Preserving Planetary-Event Null
==========================================================================

Control for the seasonal confound in the Step 2.6 planetary-event test.

MOTIVATION
==========
Step 2.6 detects elevated Gaussian-pulse significance in ±120-day windows
centered on planetary conjunction/opposition dates (~60% of events vs ~20%
under the permutation null). The Step 2.6 null, however, shuffles daily
coherence values across dates, which destroys the seasonal autocorrelation
documented in Steps 2.3/2.8 while the real events retain their calendar
positions on the seasonal curve. A season-destroying null cannot distinguish
an event-locked anomaly from ordinary sampling of the seasonal modulation.

This step rebuilds the daily coherence series from the raw SPP NPZ products
using the identical pair-coherence methodology as Step 2.6, then applies
three season-preserving controls:

1. UNIFORM-DATE NULL (season preserved): ensembles of fake events drawn at
   uniformly random dates on the real daily series. The null now samples the
   true seasonal autocorrelation; the real event dates must exceed it to
   claim an event-locked anomaly.
2. CALENDAR-TWIN CONTROL (season matched): each real event is re-evaluated
   at the same day-of-year in the other years of the data span. Seasonal
   phase is held fixed while the planetary-event lock is broken.
3. DESEASONALIZED ANALYSIS (season removed): annual and semi-annual
   harmonics are fitted to the daily series and subtracted; the real events
   and the uniform-date null are re-scored on the residual series.

Additionally the tidal-gradient scaling channel (GM/r^3) flagged as untested
in Step 2.6 is evaluated on the stored per-event distances.

Outputs:
    - results/outputs/step_2_6b_daily_coherence_baseline_clock_bias.csv
    - results/outputs/step_2_6b_seasonal_null.json

Author: Matthew Lukin Smawfield
Date: 23 September 2026
License: CC-BY-4.0
"""

import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime, timedelta
from scipy import stats
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import os

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SCRIPT_DIR))

from scripts.utils.logger import print_status

import step_2_6_planetary_events as s26

PROCESSED_DIR = ROOT / "data" / "processed"
OUTPUTS_DIR = ROOT / "results" / "outputs"
OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)

DAILY_CACHE = OUTPUTS_DIR / "step_2_6b_daily_coherence_baseline_clock_bias.csv"
OUTPUT_JSON = OUTPUTS_DIR / "step_2_6b_seasonal_null.json"

EVENT_WINDOW_DAYS = s26.EVENT_WINDOW_DAYS
N_NULL_SETS = 200          # uniform-date null ensembles
RNG_SEED = 20260923


# ============================================================================
# PHASE 1: rebuild the daily coherence series (identical methodology to 2.6)
# ============================================================================
#
# Performance note: Step 2.6 issues one worker task per (pair, date), and each
# task rescans both stations' full DatetimeIndex (~3 years of 5-min epochs) to
# select a single day. With ~2,551 pairs x 1,096 dates that is ~2.8M full-index
# scans. Here the identical day selection is pre-grouped once per station
# (the index is time-sorted, so each day is a contiguous slice) and each worker
# task evaluates all days for one pair. The inner coherence computation is
# line-for-line identical to s26._compute_pair_coherence_for_date.

_w_station_days = {}


def _init_worker_days(payload):
    global _w_station_days
    _w_station_days = payload


def _build_day_table(series):
    """Group a time-sorted Series into contiguous per-date slices."""
    vals = series.values.astype(float)
    # int64 proleptic-ordinal day ids, equivalent to index.normalize()
    day_ids = series.index.values.astype('datetime64[D]').astype(np.int64)
    uniq, starts = np.unique(day_ids, return_index=True)
    ends = np.append(starts[1:], len(day_ids))
    return {'vals': vals, 'starts': starts, 'ends': ends,
            'ord2i': {int(o): i for i, o in enumerate(uniq)}}


def _pair_day_core(day_ts1, day_ts2):
    """Inner per-day MSC + phase-alignment computation.

    Line-for-line identical to the body of
    s26._compute_pair_coherence_for_date (lines following the day mask).
    """
    min_len = min(len(day_ts1), len(day_ts2))
    if min_len < 50:
        return None

    day_ts1 = day_ts1[:min_len]
    day_ts2 = day_ts2[:min_len]

    valid_mask = np.isfinite(day_ts1) & np.isfinite(day_ts2)
    if np.sum(valid_mask) < 50:
        return None

    s1 = day_ts1[valid_mask]
    s2 = day_ts2[valid_mask]

    t = np.arange(len(s1))
    s1 = s1 - np.polyval(np.polyfit(t, s1, 1), t)
    s2 = s2 - np.polyval(np.polyfit(t, s2, 1), t)

    from scipy.signal import csd, welch
    nperseg = min(512, len(s1) // 2)
    if nperseg < 64:
        return None

    try:
        freqs, Pxy = csd(s1, s2, fs=s26.FS_HZ, nperseg=nperseg,
                         detrend='constant')
        _, Pxx = welch(s1, fs=s26.FS_HZ, nperseg=nperseg, detrend='constant')
        _, Pyy = welch(s2, fs=s26.FS_HZ, nperseg=nperseg, detrend='constant')

        band_mask = (freqs > 0) & (freqs >= s26.F1_HZ) & (freqs <= s26.F2_HZ)
        if not np.any(band_mask):
            return None

        Pxy_band = Pxy[band_mask]
        denom = Pxx[band_mask] * Pyy[band_mask]
        valid = denom > 0
        if not np.any(valid):
            return None

        coh_squared = np.abs(Pxy_band[valid]) ** 2 / denom[valid]
        magnitudes = np.sqrt(coh_squared)
        mean_coherence = float(np.mean(magnitudes))

        phases = np.angle(Pxy_band[valid])
        if len(magnitudes) > 0 and np.sum(magnitudes) > 0:
            complex_phases = np.exp(1j * phases)
            weighted_complex = np.average(complex_phases, weights=magnitudes)
            weighted_phase = np.angle(weighted_complex)
            phase_alignment = float(np.cos(weighted_phase))
        else:
            phase_alignment = np.nan

        return {'coherence': mean_coherence, 'phase_alignment': phase_alignment}
    except Exception:
        return None


def _compute_pair_all_days(args):
    """Compute per-day coherence for one pair across all requested dates.

    Returns a list of records identical in fields to
    s26._compute_pair_coherence_for_date output.
    """
    sta1, sta2, dist, date_ords = args
    A = _w_station_days[sta1]
    B = _w_station_days[sta2]
    epoch = datetime(1970, 1, 1).toordinal()
    out = []
    for do in date_ords:
        i = A['ord2i'].get(do)
        j = B['ord2i'].get(do)
        if i is None or j is None:
            continue
        r = _pair_day_core(A['vals'][A['starts'][i]:A['ends'][i]],
                           B['vals'][B['starts'][j]:B['ends'][j]])
        if r is None:
            continue
        dt = datetime.fromordinal(int(do) + epoch)
        out.append({'year': dt.year, 'doy': dt.timetuple().tm_yday,
                    'date': pd.Timestamp(dt), 'dist': dist,
                    'coherence': r['coherence'],
                    'phase_alignment': r['phase_alignment']})
    return out


def build_daily_series():
    """Rebuild the daily mean coherence series from raw NPZ products."""
    if DAILY_CACHE.exists():
        print_status(f"  Loading cached daily series: {DAILY_CACHE.name}", "INFO")
        df = pd.read_csv(DAILY_CACHE, parse_dates=['date'])
        return df

    with open(PROCESSED_DIR / "station_coordinates.json") as f:
        coords_ecef = json.load(f)
    station_coords = {}
    for sta, ecef in coords_ecef.items():
        lat, lon = s26.ecef_to_lla(ecef[0], ecef[1], ecef[2])
        station_coords[sta] = {'lat': lat, 'lon': lon}

    station_filter = os.environ.get('STATION_FILTER', 'optimal_100_metadata.json')
    with open(PROCESSED_DIR / station_filter) as f:
        good = set(json.load(f)['stations'])
    print_status(f"  Station filter: {station_filter} ({len(good)} stations)", "INFO")

    npz_files = sorted(PROCESSED_DIR.glob("*.npz"))
    station_files, available_dates = {}, set()
    for f in npz_files:
        parts = f.stem.split('_')
        if len(parts) < 2:
            continue
        sta = parts[0]
        if sta not in good or sta not in station_coords:
            continue
        station_files.setdefault(sta, []).append(f)
        try:
            ds = parts[1]
            if len(ds) == 7:
                y, doy = int(ds[:4]), int(ds[4:])
                available_dates.add(pd.Timestamp(datetime(y, 1, 1) + timedelta(days=doy - 1)))
        except Exception:
            pass

    print_status(f"  Loading {len(station_files)} stations...", "PROCESS")
    from scripts.utils.data_alignment import load_aligned_data
    station_data = {}
    for i, (sta, files) in enumerate(station_files.items()):
        if i % 10 == 0:
            print(f"\r    {i}/{len(station_files)}", end="", flush=True)
        try:
            out = load_aligned_data(files, year=None, keys=['clock_bias_ns'])
            if isinstance(out, pd.DataFrame):
                if 'clock_bias_ns' not in out.columns:
                    continue
                out = out['clock_bias_ns']
            series = out.dropna()
            if len(series) > 500:
                station_data[sta] = {
                    'clock_bias': series,
                    'lat': station_coords[sta]['lat'],
                    'lon': station_coords[sta]['lon'],
                }
        except Exception:
            continue
    print()
    print_status(f"  Loaded {len(station_data)} stations", "SUCCESS")

    # Pair list (identical distance cuts to Step 2.6)
    stations = list(station_data.keys())
    pair_list = []
    for i, s1 in enumerate(stations):
        for s2 in stations[i + 1:]:
            d = s26.haversine(station_data[s1]['lat'], station_data[s1]['lon'],
                              station_data[s2]['lat'], station_data[s2]['lon'])
            if s26.MIN_DISTANCE_KM <= d <= s26.MAX_DISTANCE_KM:
                pair_list.append((s1, s2, d))
    dates = sorted(available_dates)
    print_status(f"  {len(pair_list)} pairs x {len(dates)} dates", "INFO")

    # Pre-group each station's series into per-date slices (one-time O(N) pass),
    # then evaluate one task per pair instead of per (pair, date).
    station_days = {sta: _build_day_table(d['clock_bias'])
                    for sta, d in station_data.items()}
    date_ords = np.array(
        [d.date().toordinal() - datetime(1970, 1, 1).toordinal() for d in dates],
        dtype=np.int64)

    tasks = [(s1, s2, d, date_ords) for s1, s2, d in pair_list]
    n_workers = max(1, multiprocessing.cpu_count() - 1)
    print_status(f"  Computing pair-day coherences: {len(tasks)} pairs x "
                 f"{len(date_ords)} days ({n_workers} workers)...", "PROCESS")

    results = []
    with ProcessPoolExecutor(max_workers=n_workers,
                             initializer=_init_worker_days,
                             initargs=(station_days,)) as ex:
        futs = [ex.submit(_compute_pair_all_days, t) for t in tasks]
        done = 0
        for fut in as_completed(futs):
            results.extend(fut.result())
            done += 1
            if done % 100 == 0 or done == len(tasks):
                print(f"\r    {done}/{len(tasks)} pairs", end="", flush=True)
    print()
    print_status(f"  {len(results)} valid pair-day coherences", "SUCCESS")

    agg_msc = defaultdict(list)
    agg_pa = defaultdict(list)
    for r in results:
        key = (r['year'], r['doy'])
        agg_msc[key].append(r['coherence'])
        if np.isfinite(r['phase_alignment']):
            agg_pa[key].append(r['phase_alignment'])

    rows = []
    for (y, doy), vals in sorted(agg_msc.items()):
        if len(vals) < s26.MIN_DAILY_PAIRS:
            continue
        rows.append({
            'year': y, 'doy': doy,
            'date': pd.Timestamp(datetime(y, 1, 1) + timedelta(days=doy - 1)),
            'mean_coherence': float(np.mean(vals)),
            'std_coherence': float(np.std(vals)),
            'n_pairs': len(vals),
            'phase_alignment': float(np.mean(agg_pa.get((y, doy), [np.nan]))),
        })
    daily_df = pd.DataFrame(rows).sort_values('date').reset_index(drop=True)
    daily_df.to_csv(DAILY_CACHE, index=False)
    print_status(f"  Daily series: {len(daily_df)} days -> {DAILY_CACHE.name}", "SUCCESS")
    return daily_df


# ============================================================================
# PHASE 2: season-preserving null tests
# ============================================================================

def score_event(daily_df, metric_col, year, doy):
    """Score one (year, doy) through the Step 2.6 Gaussian-pulse pipeline."""
    df = daily_df[['date', metric_col]].rename(columns={metric_col: 'mean_coherence'})
    res = s26.analyze_planetary_event(df, year, doy, EVENT_WINDOW_DAYS)
    if not res['success']:
        return None
    if res.get('gaussian_fit'):
        return {'sigma_level': res['gaussian_fit']['sigma_level'],
                'is_significant': res['gaussian_fit']['is_significant']}
    sa = res.get('simple_analysis', {})
    return {'sigma_level': sa.get('sigma_level', 0.0),
            'is_significant': sa.get('is_significant', False)}


def score_event_set(daily_df, metric_col, events):
    """Score a list of (year, doy) events; return rate, sigmas."""
    n_sig, sigmas = 0, []
    for y, doy in events:
        r = score_event(daily_df, metric_col, y, doy)
        if r is None:
            continue
        sigmas.append(r['sigma_level'])
        n_sig += r['is_significant']
    n = len(sigmas)
    return {'n': n, 'n_sig': n_sig,
            'rate': n_sig / n if n else np.nan,
            'sigmas': sigmas,
            'mean_sigma': float(np.mean(sigmas)) if sigmas else np.nan}


def uniform_date_null(daily_df, metric_col, n_events, n_sets, rng):
    """Draw fake-event ensembles at uniform random dates on the REAL series."""
    lo = daily_df['date'].min() + pd.Timedelta(days=EVENT_WINDOW_DAYS)
    hi = daily_df['date'].max() - pd.Timedelta(days=EVENT_WINDOW_DAYS)
    pool = daily_df[(daily_df['date'] >= lo) & (daily_df['date'] <= hi)]
    cand = [(int(r.year), int(r.doy)) for r in pool.itertuples()]
    rates, sigmas_all = [], []
    for _ in range(n_sets):
        idx = rng.choice(len(cand), size=n_events, replace=False)
        ev = [cand[i] for i in idx]
        sc = score_event_set(daily_df, metric_col, ev)
        rates.append(sc['rate'])
        sigmas_all.extend(sc['sigmas'])
    return {'rates': np.array(rates), 'sigmas': np.array(sigmas_all),
            'n_candidates': len(cand)}


def deseasonalize(daily_df, metric_col, n_harmonics=2):
    """Remove annual + semi-annual harmonics; return residual DataFrame."""
    t = (daily_df['date'] - daily_df['date'].min()).dt.days.values.astype(float)
    y = daily_df[metric_col].values
    cols = [np.ones_like(t)]
    for k in range(1, n_harmonics + 1):
        cols += [np.sin(2 * np.pi * k * t / 365.25), np.cos(2 * np.pi * k * t / 365.25)]
    X = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    # lstsq's LAPACK solve can leave sticky FP flags that flush as spurious
    # RuntimeWarnings on the next ufunc; suppress and clear them here.
    with np.errstate(all='ignore'):
        resid = y - X @ coef
    out = daily_df.copy()
    out[metric_col] = resid
    return out, coef


def get_real_events():
    events = []
    for year, planets in s26.PLANETARY_EVENTS_BY_YEAR.items():
        for planet, ev in planets.items():
            for etype, doys in ev.items():
                if etype != 'type':
                    for doy in doys:
                        events.append({'planet': planet, 'year': year, 'doy': doy})
    return events


def gm_r3_test(events_scored):
    """Tidal-gradient scaling: sigma_level vs GM/r^3 (JPL distances)."""
    GM = {'Mercury': 22032, 'Venus': 324860, 'Mars': 42828,
          'Jupiter': 126690000, 'Saturn': 37931000}
    xs, ys = [], []
    for e in events_scored:
        if e.get('sigma_level') is None or e.get('distance_au') in (None, 0):
            continue
        xs.append(GM[e['planet']] / e['distance_au'] ** 3)
        ys.append(e['sigma_level'])
    if len(xs) < 5:
        return {'valid': False}
    r, p = stats.pearsonr(np.log10(xs), ys)
    return {'valid': True, 'n': len(xs), 'pearson_r': float(r), 'p_value': float(p)}


def stored_event_distances():
    """Pull per-event distance_au from the Step 2.6 primary output."""
    f = OUTPUTS_DIR / "step_2_6_planetary_events_baseline_clock_bias_msc.json"
    if not f.exists():
        return {}
    d = json.load(open(f))
    out = {}
    for planet, res in d.get('planetary_results', {}).items():
        for e in res.get('events', []):
            out[(planet, e['year'], e['doy'])] = e.get('distance_au')
    return out


def run_metric_tests(daily_df, metric_col, real_events, rng, label):
    print_status(f"  --- {label} ---", "INFO")
    real_sc = score_event_set(daily_df, metric_col,
                              [(e['year'], e['doy']) for e in real_events])
    print_status(f"    real events: {real_sc['n_sig']}/{real_sc['n']} "
                 f"({100 * real_sc['rate']:.1f}%), mean sigma {real_sc['mean_sigma']:.2f}", "INFO")

    # Uniform-date null on the real (season-preserving) series
    null = uniform_date_null(daily_df, metric_col, real_sc['n'], N_NULL_SETS, rng)
    null_rates = null['rates']
    p_rate = float(np.mean(null_rates >= real_sc['rate']))
    stat, p_mw = stats.mannwhitneyu(real_sc['sigmas'], null['sigmas'],
                                  alternative='greater')
    print_status(f"    uniform-date null: {100 * np.mean(null_rates):.1f}% "
                 f"+/- {100 * np.std(null_rates):.1f}% over {N_NULL_SETS} sets; "
                 f"p(rate)={p_rate:.3f}, MW p={p_mw:.4f}", "INFO")

    # Calendar-twin control: same DOY, other years
    years = sorted(daily_df['year'].unique())
    twin_events = [(y2, e['doy']) for e in real_events for y2 in years if y2 != e['year']]
    twin_sc = score_event_set(daily_df, metric_col, twin_events)
    print_status(f"    calendar twins (same DOY, other years): {twin_sc['n_sig']}/{twin_sc['n']} "
                 f"({100 * twin_sc['rate']:.1f}%), mean sigma {twin_sc['mean_sigma']:.2f}", "INFO")

    # Deseasonalized re-analysis
    resid_df, coef = deseasonalize(daily_df, metric_col)
    real_resid = score_event_set(resid_df, metric_col,
                                 [(e['year'], e['doy']) for e in real_events])
    null_resid = uniform_date_null(resid_df, metric_col, real_resid['n'], N_NULL_SETS, rng)
    p_rate_resid = float(np.mean(null_resid['rates'] >= real_resid['rate']))
    stat_r, p_mw_r = stats.mannwhitneyu(real_resid['sigmas'], null_resid['sigmas'],
                                      alternative='greater')
    print_status(f"    deseasonalized: real {100 * real_resid['rate']:.1f}% vs null "
                 f"{100 * np.mean(null_resid['rates']):.1f}%; "
                 f"p(rate)={p_rate_resid:.3f}, MW p={p_mw_r:.4f}", "INFO")

    return {
        'metric': metric_col,
        'real': {'n': real_sc['n'], 'n_sig': real_sc['n_sig'],
                 'rate': real_sc['rate'], 'mean_sigma': real_sc['mean_sigma'],
                 'sigmas': real_sc['sigmas']},
        'uniform_date_null': {
            'n_sets': N_NULL_SETS,
            'mean_rate': float(np.mean(null_rates)),
            'std_rate': float(np.std(null_rates)),
            'p_rate': p_rate,
            'mannwhitney_p': float(p_mw),
            'mean_sigma': float(np.mean(null['sigmas'])),
        },
        'calendar_twins': {'n': twin_sc['n'], 'n_sig': twin_sc['n_sig'],
                           'rate': twin_sc['rate'], 'mean_sigma': twin_sc['mean_sigma']},
        'deseasonalized': {
            'harmonic_coeffs': [float(c) for c in coef],
            'real': {'n': real_resid['n'], 'n_sig': real_resid['n_sig'],
                     'rate': real_resid['rate'], 'mean_sigma': real_resid['mean_sigma']},
            'null': {'mean_rate': float(np.mean(null_resid['rates'])),
                     'std_rate': float(np.std(null_resid['rates'])),
                     'p_rate': p_rate_resid,
                     'mannwhitney_p': float(p_mw_r)},
        },
    }


def main():
    print_status("=" * 80, "INFO")
    print_status("STEP 2.6b: Season-Preserving Planetary-Event Null", "INFO")
    print_status("=" * 80, "INFO")

    rng = np.random.default_rng(RNG_SEED)

    print_status("[1/3] Building daily coherence series (baseline clock_bias)...", "PROCESS")
    daily_df = build_daily_series()
    span = (daily_df['date'].min(), daily_df['date'].max())
    print_status(f"  Series: {len(daily_df)} days, {span[0].date()} -> {span[1].date()}, "
                 f"mean {daily_df['mean_coherence'].mean():.4f} "
                 f"+/- {daily_df['mean_coherence'].std():.4f}", "INFO")

    real_events = get_real_events()
    dist_map = stored_event_distances()

    # Seasonal placement diagnostic: coherence level at real event dates
    ev_idx = [(e['year'], e['doy']) for e in real_events]
    lookup = {(int(r.year), int(r.doy)): r.mean_coherence for r in daily_df.itertuples()}
    ev_coh = np.array([lookup[k] for k in ev_idx if k in lookup])
    z = (ev_coh.mean() - daily_df['mean_coherence'].mean()) / daily_df['mean_coherence'].std()
    print_status(f"  Event-date coherence placement: z = {z:+.2f} "
                 f"(mean at events {ev_coh.mean():.4f} vs series {daily_df['mean_coherence'].mean():.4f})", "INFO")

    print_status("[2/3] Season-preserving nulls...", "PROCESS")
    results = {}
    for metric_col in ['mean_coherence', 'phase_alignment']:
        results[metric_col] = run_metric_tests(daily_df, metric_col, real_events,
                                             rng, f"baseline clock_bias / {metric_col}")

    print_status("[3/3] Tidal-gradient (GM/r^3) scaling on stored event data...", "PROCESS")
    stored = json.load(open(OUTPUTS_DIR / "step_2_6_planetary_events_baseline_clock_bias_msc.json"))
    scored = []
    for planet, res in stored.get('planetary_results', {}).items():
        for e in res.get('events', []):
            scored.append({'planet': planet, 'sigma_level': e.get('sigma_level'),
                           'distance_au': e.get('distance_au'),
                           'modulation_pct': e.get('modulation_pct')})
    gmr3_sigma = gm_r3_test(scored)
    gmr3_mod = gm_r3_test([{**s, 'sigma_level': abs(s['modulation_pct'] or 0)} for s in scored])
    print_status(f"    sigma vs GM/r^3: r={gmr3_sigma.get('pearson_r', float('nan')):+.3f}, "
                 f"p={gmr3_sigma.get('p_value', float('nan')):.4f}", "INFO")
    print_status(f"    |mod| vs GM/r^3: r={gmr3_mod.get('pearson_r', float('nan')):+.3f}, "
                 f"p={gmr3_mod.get('p_value', float('nan')):.4f}", "INFO")

    output = {
        'step': '2.6b',
        'name': 'season_preserving_planetary_event_null',
        'timestamp': datetime.now().isoformat(),
        'methodology': ('Rebuilds the Step 2.6 daily coherence series from raw NPZ, then '
                        'applies season-preserving controls: uniform-date null on the real '
                        'series, same-DOY other-year calendar twins, and harmonic-'
                        'deseasonalized rescoring; plus GM/r^3 tidal-gradient scaling.'),
        'daily_series': {
            'n_days': len(daily_df),
            'date_min': str(span[0].date()), 'date_max': str(span[1].date()),
            'mean_coherence': float(daily_df['mean_coherence'].mean()),
            'std_coherence': float(daily_df['mean_coherence'].std()),
            'event_date_z': float(z),
            'cache': DAILY_CACHE.name,
        },
        'event_window_days': EVENT_WINDOW_DAYS,
        'n_null_sets': N_NULL_SETS,
        'tests': results,
        'gm_r3_scaling': {'sigma_vs_gm_r3': gmr3_sigma, 'mod_vs_gm_r3': gmr3_mod},
    }
    # strip sigma lists for JSON compactness
    for m in output['tests'].values():
        m['real'].pop('sigmas', None)

    with open(OUTPUT_JSON, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    print_status(f"Results saved: {OUTPUT_JSON}", "SUCCESS")


if __name__ == "__main__":
    main()
