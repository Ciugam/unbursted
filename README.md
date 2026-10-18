# Unbursted

Recovering the true number of free molecules in FCS measurements with free and
DNA-bound EtBr, using a threshold method that rejects the bursts produced by
the bright intercalated species. What remains — the "unbursted" trace — is
analysed via autocorrelation to recover N_free and the diffusion time.

> Note: all parameters reported below refer to the project *as currently set
> up*. Every value can be varied if desired (diffusion coefficients, N_avg,
> bound fractions, box/PSF geometry, timing, etc.).

## Objective

FCS measurements are polluted by bursts of high fluorescence produced by bright
species (e.g. EtBr intercalated in DNA). This project aims to:

- build a simulation pipeline (Brownian motion → timestamps) reproducing free
  and DNA-bound EtBr signals;
- implement a **threshold method** that rejects burst windows and analyses the
  remaining trace via autocorrelation (1-species model) to recover the correct
  number of **free** molecules and the diffusion time.

Free EtBr is the reference case (same bright fluorophore as the DNA-bound one,
used to calibrate the free-species behaviour).

## Pipeline (3 phases)

### Phase 1 — Brownian motion simulation
Notebook `1_SIM_brownian_motion.ipynb` (done). Uses `pybromo` to simulate
particle diffusion in a confocal volume and produces fluorescence emission
traces saved to `Pybromo_simulations/` (*.hdf5).

- Box: 5 × 5 × 10 µm (reflective boundaries)
- PSF: 3D Gaussian, w0 = 0.25 µm, z0 = 1.0 µm (ellipsoid confocal volume)
- Free EtBr: D = 425 µm²/s
- DNA-bound EtBr modelled as whole DNA blobs (coils): D = 1 µm²/s
- n_avg = 1 free molecule in the confocal volume; blob populations
  simulated: 0 (free only), 0.12 and 0.18 blobs in the confocal volume
- Absolute particle count in the box derived from volume/concentration:
  N_box = n_avg · (V_box / V_conf); for blobs, at least one is always
  generated (N_blobs = max(1, int(n_avg_bound · V_box / V_conf)))
- Timing: t_step = 5 µs, t_max = 60 s (all runs)

### Phase 2 — Timestamp generation
Notebook `2_SIM_timestamps_gen.ipynb` (done). Converts Phase 1 emission
traces into photon arrival timestamps, saved to `Pybromo_timestamps/` (*.hdf5).

- Emission rates (identical across all runs): max_rate_free = 600 cps,
  max_rate_dna = 24000 cps, background = 275 + 2400 cps

### Phase 3 — Data analysis
Notebook `3_ANALYSIS.ipynb` (**done**).
- Autocorrelation of simulated timestamps, 1-species and double-population fits, compare N and diffusion time vs input values
- Cross-correlation of real data (detectors m1 & m2), 1-species fit
- Threshold method implemented: median + k·MAD threshold with margin_bins sweep
- Centralized diffusion/fit functions in `toolkit.py`: diffusion_3d_single/double/single_offset, fit_G
- Validation sections for Free EtBr reference, EtBr-DNA 0.06 blobs Sol03 vs ac470b, EtBr-DNA 0.08 blobs Sol05 vs 00512f
- DataFrame summaries for Real vs Sim with τ_diff, D, A0_free, N_eff
- k-sweep plots and photons remaining vs margin_bins for k=0.5/1.0/1.5
- Unified English summary with tables and polydispersity note
- Sweet spot chosen: k=1.0, margin=0 for both DNA concentrations

## Repository structure

    unbursted/
    ├── 1_SIM_brownian_motion.ipynb    # Phase 1: Brownian motion simulation ✓
    ├── 2_SIM_timestamps_gen.ipynb     # Phase 2: timestamp generation ✓
    ├── 3_ANALYSIS.ipynb               # Phase 3: data analysis ✓
    ├── toolkit.py                     # Centralized diffusion/fit functions
    ├── Pybromo_simulations/           # Phase 1 output (emission traces, *.hdf5)
    ├── Pybromo_timestamps/            # Phase 2 output (timestamps, *.hdf5)
    ├── Real_data/                     # Experimental data (.asc)
    │   ├── m1/                        # Detector 1 (same runs as m2)
    │   └── m2/                        # Detector 2 (same runs as m1)
    ├── docs/
    │   └── EtBr_DNA_Samples.xlsx
    ├── .gitignore
    └── README.md

Real data: files are named `SolXX_mY.asc` (run Sol01, 03, 05, 14; detector
m1/m2); both detectors acquire simultaneously in every run.

Simulation files embed the main specs in the name, e.g.
`pybromo_…_P954_D4.2e-10_P114_D1e-12_7093pM_step5.0us_t_max60.0s_ID0-0.hdf5`
(P/D = particle count and diffusion coefficient in m²/s per population,
concentration in pM, step/t_max = time binning, ID = chunk index).

## Completed items

- Timestamp histogram binning aligned with the experimental one for all runs
- Threshold method validated on Free EtBr and EtBr-DNA 0.06/0.08 blob simulations
- Sweet spot k=1.0, margin=0 adopted for both DNA concentrations

## Author

- **Ciugam** (theroot@ciugam.net)
