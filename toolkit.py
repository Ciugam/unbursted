import tables
import os
import numpy as np 
import pycorrelate as pyc
from scipy.optimize import curve_fit

def load_sim_timestamps(path_ts, fname, return_par=False):
    """Load photon timestamps from a PyTables simulation file.

    Args:
        path_ts: directory containing the .hdf5 file.
        fname: name of the .hdf5 file.
        return_par: if True, also return the particle-ID array (the
            background photons carry the ID `num_particles`, i.e. last
            index + 1). Useful to separate signal from background before
            computing the ACF.

    Returns:
        numpy array of timestamps in seconds (and, if return_par, the
        particle-ID array).
    """

    with tables.open_file(os.path.join(path_ts, fname), mode='r') as f:
        name = [n for n in f.root.timestamps._v_children
                if not n.endswith('_par')][0]     # first timestamps dataset
        ts = f.root.timestamps[name]
        ts_sim = ts.read()
        # simulated clock resolution picked up directly from dataset
        clk_p = ts.attrs['clk_p'] 
        ts_sim = ts_sim * clk_p 
        if return_par:
            par = f.root.timestamps[name + '_par'].read()
            return ts_sim, par
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
        tau, G, Y where G is normalized correlation and Y is the number of
        photon pairs in each lag bin (so sigma = 1/sqrt(Y) is the correct
        Poisson estimate of the uncertainty on G).
    """
    t1 = np.round(np.sort(ts1) / unit).astype(np.int64)
    t2 = np.round(np.sort(ts2) / unit).astype(np.int64)
    exp_max = int(np.ceil(np.log10(t_max / unit)))
    bins = pyc.make_loglags(0, exp_max, bins_per_dec)
    G = pyc.pcorrelate(t1, t2, bins=bins, normalize=True)
    Y_raw = pyc.pcorrelate(t1, t2, bins=bins, normalize=False)  # pairs / bin_width (ticks)
    tau = 0.5 * (bins[1:] + bins[:-1]) * unit

    # Numero effettivo di coppie di fotoni per bin di lag: con questa Y,
    # sigma = 1/sqrt(Y) e' la stima di Poisson corretta dell'incertezza di G.
    dbin = np.diff(bins).astype(np.float64)
    Y = Y_raw * dbin
    return tau, G, Y

def unburst(ts_array, bins, threshold):
    """Remove burst windows from a timestamp array.

    Args:
        ts_array: 1D array of timestamps.
        bins: bin edges used to detect bursts.
        threshold: count threshold for burst detection.

    Returns:
        filtered timestamp array with burst bins removed.
    """
    counts = np.histogram(ts_array, bins)[0]
    above = counts > threshold
    padded = np.concatenate(([False], above, [False]))
    diffs = np.diff(padded.astype(int))
    starts, ends = np.where(diffs == 1)[0], np.where(diffs == -1)[0]
    flag = np.zeros_like(above, bool)
    for s, e in zip(starts, ends):
        flag[s - 1:e] = True
    bins_digitized = np.digitize(ts_array, bins)
    unburst_mask = ~np.isin(bins_digitized - 1, np.where(flag)[0])
    return ts_array[unburst_mask]

# Default confocal PSF (W0, Z0) in meters. These are used only by the
# convenience wrappers below; fits should pass W0/Z0 explicitly (e.g. read
# from the simulation metadata).
W0_DEFAULT = 0.25e-6
Z0_DEFAULT = 1.0e-6


def _G3d(tau, tau_diff, A0, W0, Z0):
    KAPPA = Z0 / W0
    return A0 / (1 + tau/tau_diff) / np.sqrt(1 + tau/(KAPPA**2 * tau_diff))


def diffusion_3d_single(tau, tau_diff, A0, W0=W0_DEFAULT, Z0=Z0_DEFAULT):
    """Single-species 3D-diffusion ACF: G = 1 + A0 * f(tau)."""
    return 1 + _G3d(tau, tau_diff, A0, W0, Z0)


def diffusion_3d_double(tau, tau_diff_1, A0_1, tau_diff_2, A0_2, W0=W0_DEFAULT, Z0=Z0_DEFAULT):
    """Two-species 3D-diffusion ACF: G = 1 + A0_1*f1 + A0_2*f2."""
    return 1 + _G3d(tau, tau_diff_1, A0_1, W0, Z0) + _G3d(tau, tau_diff_2, A0_2, W0, Z0)





def fit_G(tau, G, guesses, bounds, model="single", tau_range=(1e-5, 2e-4), label="",
          sigma=None, W0=W0_DEFAULT, Z0=Z0_DEFAULT, kappa_fixed=True, verbose=False):
    """Fit the ACF on a restricted tau range.

    Args:
        tau_range: tuple (tau_lo, tau_hi) in seconds. Only points with
            tau_lo <= tau <= tau_hi are used for the fit. Choose it from the
            ACF plot: below ~5-10 us the curve is noisy / affected by short-
            time artifacts, above ~10*tau_diff the decay is dominated by the
            long-time artifact (G falling artificially).
        W0, Z0: confocal PSF radii in meters, FIXED during the fit. For
            simulation data read them from the file metadata; for real data
            pass your calibrated values (or estimate W0 from a bead measurement).
        kappa_fixed: if True (default) Z0 is used as given and kappa = Z0/W0.
            If False, Z0 is fitted together with the model parameters (W0
            stays fixed). Note: W0 and Z0 are strongly degenerate in FCS;
            prefer fixing both from independent calibration.
        sigma: per-point standard deviation of G (absolute_sigma=True).

    Returns:
        (popt, pcov, W0_fit, Z0_fit). popt/pcov refer to the model parameters
        (tau_diff, A0, ...); W0_fit/Z0_fit echo the PSF values used, so
        callers can compute D = W0^2/(4*tau_diff) with the correct W0.
    """
    n_par = {"single": 2, "double": 4}
    if model not in n_par:
        raise ValueError(f'model "{model}" non valida: usa "single" o "double"')
    n_par = n_par[model]

    # Build the fitting function with W0/Z0 as extra trailing parameters
    def _wrap(tau_, *p):
        pars_model = p[:n_par]
        W0_ = W0  # W0 is always fixed: it is set via the W0 argument
        Z0_ = Z0 if kappa_fixed else p[n_par]
        if model == "single":
            return diffusion_3d_single(tau_, *pars_model, W0=W0_, Z0=Z0_)
        else:
            return diffusion_3d_double(tau_, *pars_model, W0=W0_, Z0=Z0_)

    # Extend guesses and bounds with the PSF parameters.
    # W0 is fixed unless bounds contain an extra pair (len > n_par):
    # in that case the last pair is used as (W0_lo, W0_hi) and W0 is fitted.
    p0 = list(guesses)
    lo, hi = list(bounds[0]), list(bounds[1])
    if not kappa_fixed:
        p0 += [Z0]; lo += [Z0 * 0.5]; hi += [Z0 * 2.0]

    mask = (tau >= tau_range[0]) & (tau <= tau_range[1])
    sigma_m = sigma[mask] if sigma is not None else None
    popt_all, pcov_all = curve_fit(_wrap, tau[mask], G[mask], p0=p0,
                                   bounds=(lo, hi), sigma=sigma_m, absolute_sigma=True)
    popt = popt_all[:n_par]
    pcov = pcov_all[:n_par, :n_par]
    W0_fit = W0
    Z0_fit = Z0 if kappa_fixed else popt_all[n_par]

    if verbose:
        names = {"single": ["tau_diff", "A0"],
                 "double": ["tau1", "A0_1", "tau2", "A0_2"]}[model]
        perr = np.sqrt(np.diag(pcov))
        vals = "   ".join(f"{n} = {v:.4g} +/- {e:.3g}" for n, v, e in zip(names, popt, perr))
        print(f"[{label}] {vals}   W0 = {W0_fit:.3g} m")

    return popt, pcov, W0_fit, Z0_fit

def fit_k_sweep(N_bound_vals, k_vals):
    """Fit the hyperbolic relation k = a/N + b between optimal threshold and
    number of bound molecules in the confocal volume.

    Parameters
    ----------
    N_bound_vals : array-like
        Number of bound molecules (in confocal volume) for each simulation.
    k_vals : array-like
        Corresponding optimal k values from the sweep.

    Returns
    -------
    popt : tuple (a, b)
        Fitted parameters: k = a/N + b
    pcov : ndarray
        Covariance matrix of the fit.
    """
    def hyperbolic(N, a, b):
        return a / N + b

    popt, pcov = curve_fit(hyperbolic, N_bound_vals, k_vals, p0=(1.0, -0.5))
    return popt, pcov
