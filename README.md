# Unbursted

Recovering the true number of free molecules in FCS measurements with free and
DNA-bound EtBr, using a threshold method that rejects the bursts produced by
the bright intercalated species. What remains — the "unbursted" trace — is
analysed via autocorrelation to recover N_free and the diffusion time.

## Objective

FCS measurements are polluted by bursts of high fluorescence produced by bright
species (e.g. EtBr intercalated in DNA). This project aims to:

- build a simulation pipeline (Brownian motion → timestamps) reproducing free
  and DNA-bound EtBr signals;
- implement a **threshold method** that rejects burst windows and analyses the
  remaining trace via autocorrelation (1-species model) to recover the correct
  number of **free** molecules and the diffusion time;
- validate the workflow on noise-free simulations, then test it on real
  experimental data.

Free EtBr is the reference case (same bright fluorophore as the DNA-bound one,
used to calibrate the free-species behaviour).

## Pipeline (3 phases)

### Phase 1 — Brownian motion simulation
Notebook `1_SIM_brownian_motion.ipynb`. Uses `pybromo` to simulate particle
diffusion in a confocal volume and produces fluorescence emission traces saved
to `Pybromo_simulations/` (*.hdf5). All parameters are read from
`parameters.toml`.

- Box: 2 × 2 × 4 µm (reflective boundaries)
- PSF: 3D Gaussian, w0 = 0.30 µm, z0 = 1.20 µm (ellipsoid confocal volume)
- Free EtBr: D = 425 µm²/s, N = 100 particles in the box
- DNA-bound EtBr modelled as whole DNA blobs (coils): D = 0.5 µm²/s; blob
  populations simulated: 6, 8, 10, 12 and 15 complexes (plus a free-EtBr-only
  run with no bound species)
- Timing: t_step = 1 µs, t_max = 120 s (all runs)

### Phase 2 — Timestamp generation
Notebook `2_SIM_timestamps_gen.ipynb`. Converts Phase 1 emission traces into
photon arrival timestamps, saved to `Pybromo_timestamps/` (*.hdf5).

- Emission rates: max_rate_free = 600 cps; max_rate_dna = 10 000–12 000 cps
  (depending on blob count); background = 275 cps

### Phase 3 — Data analysis
Two dedicated notebooks:

- `3a_simulation_analysis.ipynb` — unburst workflow on the noise-free
  simulations: autocorrelation, k-sweep per blob density, N_free recovery and
  scaling law for the optimal threshold.
- `3b_realdata_analysis.ipynb` — workflow on real PicoQuant data: free-EtBr
  reference fit, then unburst vs raw-data analysis on EtBr + EtBr-DNA samples.

Shared diffusion and fitting functions live in `toolkit.py`.

## Repository structure

    unbursted/
    ├── parameters.toml                  # Simulation parameters (single source of truth)
    ├── 1_SIM_brownian_motion.ipynb      # Phase 1: Brownian motion simulation
    ├── 2_SIM_timestamps_gen.ipynb       # Phase 2: timestamp generation
    ├── 3a_simulation_analysis.ipynb     # Phase 3a: analysis of simulations
    ├── 3b_realdata_analysis.ipynb       # Phase 3b: analysis of real data
    ├── toolkit.py                       # Centralized diffusion/fit functions
    ├── Pybromo_simulations/             # Phase 1 output (emission traces, *.hdf5)
    ├── Pybromo_timestamps/              # Phase 2 output (timestamps, *.hdf5)
    ├── Real_data/                       # Experimental data (.asc)
    │   ├── m1/                          # Detector 1 (same runs as m2)
    │   └── m2/                          # Detector 2 (same runs as m1)
    ├── docs/
    │   └── EtBr_DNA_Samples.xlsx
    ├── .gitignore
    └── README.md

Real data: files are named `SolXX_mY.asc` (run Sol01, 03, 05, 14; detector
m1/m2); both detectors acquire simultaneously in every run.

Simulation files embed the main specs in the name, e.g.
`times_767b38_P100_D4.2e-10_P6_D5e-13_11001pM_step1.0us_t_max120.0s_ID0-0.hdf5`
(P/D = particle count and diffusion coefficient in m²/s per population,
concentration in pM, step/t_max = time binning, ID = chunk index).

## Author

- **Ciugam** (theroot@ciugam.net)
