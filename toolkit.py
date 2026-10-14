import tables
import os
import numpy as np 
import pycorrelate as pyc
from scipy.optimize import curve_fit

def load_sim_timestamps(path_ts, fname):
    """Load photon timestamps from a PyTables simulation file.

    Args:
        path_ts: directory containing the .hdf5 file.
        fname: name of the .hdf5 file.

    Returns:
        numpy array of timestamps in seconds.
    """

    with tables.open_file(os.path.join(path_ts, fname), mode='r') as f:
        name = [n for n in f.root.timestamps._v_children
                if not n.endswith('_par')][0]     # first timestamps dataset
        ts = f.root.timestamps[name]
        ts_sim = ts.read()
        # simulated clock resolution picked up directly from dataset
        clk_p = ts.attrs['clk_p'] 
        ts_sim = ts_sim * clk_p 
    return ts_sim

def load_asc_timestamps(path_ts, fname, clk=1/(105*1e6), factor=2.44e-12):
    """Load photon timestamps from a PicoQuant .asc file.

    Args:
        path_ts: directory containing the .asc file.
        fname: name of the .asc file.
        clk: clock period in seconds.
        factor: conversion factor for the second column.

    Returns:
        numpy array of timestamps in seconds.
    """
    # PicoQuant resolution and conversion factor
    ts_raw = np.loadtxt(os.path.join(path_ts,fname), usecols=(0,1))
    ts_arr = (ts_raw[:,0] * clk + factor * ts_raw[:,1])
    return ts_arr

def trace_stats(ts, bins):
    """Compute basic statistics of a binned photon trace.

    Args:
        ts: 1D array of timestamps.
        bins: bin edges for histogramming.

    Returns:
        dict with hist_bars, ACR, median, MAD, max counts and number of >5MAD events.
    """
    ACR = ts.size/ts.max()
    c = np.histogram(ts, bins=bins)[0]
    med = np.median(c)
    MAD = 1.4826 * np.median(np.abs(c-med))
    is_peak = c > (med + 5 * MAD)
    n_env = np.sum(is_peak & ~np.roll(is_peak, 1))
    # if last is False, np.roll is counting correctly
    if is_peak[0]:
        # if last is True, recounts +1 more
        if is_peak[-1]:
            n_env += 1
    return dict(
        hist_bars = c,
        ACR = ACR,
        median = med,
        MAD = MAD,
        max = c.max(),
        envs_over_5_MAD = n_env
    )

# function: compute G(tau) from two timestamp arrays (ts1=ts2 for ACF)
def compute_G(ts1, ts2, unit, t_max, bins_per_dec=10):
    """Compute autocorrelation / cross-correlation G(tau) via pycorrelate.

    Args:
        ts1, ts2: timestamp arrays.
        unit: time resolution in seconds used to discretize timestamps.
        t_max: maximum lag.
        bins_per_dec: number of bins per decade.

    Returns:
        tau, G, Y where G is normalized correlation and Y is raw counts.
    """
    t1 = np.round(np.sort(ts1) / unit).astype(np.int64)
    t2 = np.round(np.sort(ts2) / unit).astype(np.int64)
    exp_max = int(np.ceil(np.log10(t_max / unit)))
    bins = pyc.make_loglags(0, exp_max, bins_per_dec)
    G = pyc.pcorrelate(t1, t2, bins=bins, normalize=True)
    Y = pyc.pcorrelate(t1, t2, bins=bins, normalize=False)
    tau = 0.5 * (bins[1:] + bins[:-1]) * unit
    return tau, G, Y

def burst_runs_mask(counts, thr, margin_bins=0):
    """Return a boolean mask for burst runs above threshold.

    Args:
        counts: histogram counts per bin.
        thr: threshold value.
        margin_bins: number of bins to extend before and after each burst run.

    Returns:
        boolean array True for bins belonging to a burst run.
    """
    above = counts > thr
    padded = np.concatenate(([False], above, [False]))
    diffs = np.diff(padded.astype(int))
    starts, ends = np.where(diffs == 1)[0], np.where(diffs == -1)[0]
    flagged = np.zeros_like(above, bool)
    for s, e in zip(starts, ends):
        s_m = max(s - margin_bins - 1, 0)
        e_m = min(e + margin_bins, len(flagged))
        flagged[s_m:e_m] = True
    return flagged

def unburst(ts_array, bins, threshold, margin_bins=0):
    """Remove burst windows from a timestamp array.

    Args:
        ts_array: 1D array of timestamps.
        bins: bin edges used to detect bursts.
        threshold: count threshold for burst detection.
        margin_bins: number of bins to extend before and after each burst run.

    Returns:
        filtered timestamp array with burst bins removed.
    """
    counts = np.histogram(ts_array, bins)[0]
    flag = burst_runs_mask(counts, threshold, margin_bins)
    bins_digitized = np.digitize(ts_array, bins)
    unburst_mask = ~np.isin(bins_digitized - 1, np.where(flag)[0])
    ts_array_unbursted = ts_array[unburst_mask]
    return ts_array_unbursted

def diffusion_3d_single(tau, tau_diff, A0):
    W0 = 0.25e-6
    Z0 = 1.0e-6
    KAPPA = Z0 / W0
    return 1 + A0 / (1 + tau/tau_diff) / np.sqrt(1 + tau/(KAPPA**2 * tau_diff))

def diffusion_3d_double(tau, tau_diff_1, A0_1, tau_diff_2, A0_2):
    W0 = 0.25e-6
    Z0 = 1.0e-6
    KAPPA = Z0 / W0
    term1 = A0_1 / (1 + tau/tau_diff_1) / np.sqrt(1 + tau/(KAPPA**2 * tau_diff_1))
    term2 = A0_2 / (1 + tau/tau_diff_2) / np.sqrt(1 + tau/(KAPPA**2 * tau_diff_2))
    return 1 + term1 + term2

def diffusion_3d_single_offset(tau, tau_diff, A0, C):
    W0 = 0.25e-6
    Z0 = 1.0e-6
    KAPPA = Z0 / W0
    term = A0 / (1 + tau/tau_diff) / np.sqrt(1 + tau/(KAPPA**2 * tau_diff))
    return C + term

def fit_G(tau, G, guesses, bounds, model="single", tau_min=0, label="", sigma=None, verbose=False):
    mask = tau > tau_min
    sigma_m = sigma[mask] if sigma is not None else None
    if model == "single":
        func = diffusion_3d_single
        popt, pcov = curve_fit(func, tau[mask], G[mask], p0=guesses, bounds=bounds, sigma=sigma_m, absolute_sigma=False)
        if verbose:
            print(f"[{label}] tau_diff = {popt[0]:.4g} s   A0 = {popt[1]:.4g}")
    elif model == "single_offset":
        func = diffusion_3d_single_offset
        popt, pcov = curve_fit(func, tau[mask], G[mask], p0=guesses, bounds=bounds, sigma=sigma_m, absolute_sigma=False)
        if verbose:
            print(f"[{label}] tau_diff = {popt[0]:.4g} s   A0 = {popt[1]:.4g}   C = {popt[2]:.4g}")
    elif model == "double":
        func = diffusion_3d_double
        popt, pcov = curve_fit(func, tau[mask], G[mask], p0=guesses, bounds=bounds, sigma=sigma_m, absolute_sigma=False)
        if verbose:
            print(f"[{label}] tau1 = {popt[0]:.4g} s  A0_1 = {popt[1]:.4g}"
                  f"   tau2 = {popt[2]:.4g} s  A0_2 = {popt[3]:.4g}")
    else:
        raise ValueError(f'model "{model}" non valida: usa "single", "single_offset" o "double"')
    return popt, pcov