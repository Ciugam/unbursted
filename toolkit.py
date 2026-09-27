import tables
import os
import numpy as np 
import pycorrelate as pyc
from scipy.optimize import curve_fit

def load_sim_timestamps(path_ts, fname):
    with tables.open_file(os.path.join(path_ts, fname), mode='r') as f:
        name = [n for n in f.root.timestamps._v_children
                if not n.endswith('_par')][0]     # first timestamps dataset
        ts = f.root.timestamps[name]
        ts_sim = ts.read()
        # simulated clock resolution picked up directly from dataset
        clk_p = ts.attrs['clk_p'] 
        ts_sim = ts_sim * clk_p 
    return ts_sim

def load_asc_timestamps(path_ts, fname, clk=1/(105*1e6), factor=2.44 * 1e-12):
    # PicoQuant resolution and conversion factor
    ts_raw = np.loadtxt(os.path.join(path_ts,fname), usecols=(0,1))
    ts_arr = (ts_raw[:,0] * clk + factor * ts_raw[:,1])
    return ts_arr

def trace_stats(ts, bins):
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
    t1 = np.round(np.sort(ts1) / unit).astype(np.int64)
    t2 = np.round(np.sort(ts2) / unit).astype(np.int64)
    exp_max = int(np.ceil(np.log10(t_max / unit)))
    bins = pyc.make_loglags(0, exp_max, bins_per_dec)
    G = pyc.pcorrelate(t1, t2, bins=bins, normalize=True)
    Y = pyc.pcorrelate(t1, t2, bins=bins, normalize=False)
    tau = 0.5 * (bins[1:] + bins[:-1]) * unit
    return tau, G, Y

def burst_runs_mask(counts, thr):
    above = counts > thr
    padded = np.concatenate(([False], above, [False]))
    diffs = np.diff(padded.astype(int))
    starts, ends = np.where(diffs == 1)[0], np.where(diffs == -1)[0]
    flagged = np.zeros_like(above, bool)
    for s, e in zip(starts, ends):
        flagged[s - 1:e] = True
    return flagged

def unburst(ts_array, bins, threshold):
    counts = np.histogram(ts_array, bins)[0]
    flag = burst_runs_mask(counts, threshold)
    bins_digitized = np.digitize(ts_array, bins)
    unburst_mask = ~np.isin(bins_digitized - 1, np.where(flag)[0])
    ts_array_unbursted = ts_array[unburst_mask]
    return ts_array_unbursted