# Porous-Medium Tumor Growth UQ with gPC-SG

This repository provides the concave initial-density experiment from **Test I(a)** of:

**Ning Jiang, Liu Liu, and Huimin Yu.**  
*Efficient stochastic asymptotic-preserving scheme for tumor growth models with uncertain parameters.*  
Journal of Computational and Applied Mathematics, 473 (2026), 116881.  
https://doi.org/10.1016/j.cam.2025.116881

The paper develops a stochastic asymptotic-preserving method based on generalized polynomial chaos and stochastic Galerkin projection (gPC-SG) for porous-medium tumor growth models. It investigates how uncertain initial conditions and growth parameters affect tumor evolution, and studies the connection with the Hele–Shaw limit.

This repository contains one complete experiment, including initial-data projection, PDE solvers, statistical reconstruction, reference data, and plotting routines. The computed density mean and standard deviation are compared with a stochastic collocation reference.

## Experiment: Test I(a)

The selected experiment considers a two-dimensional spatial domain with one random input. See Section 4.1, equations (29)–(30), and Figures 1–4 of the paper.

| Parameter | Setting |
|---|---|
| Random input | \(z \sim \mathcal{U}[-1,1]\) |
| Spatial domain | \([-2.2,2.2]^2\) |
| Spatial spacing | \(\Delta x=\Delta y=0.1\) |
| Time step | \(\Delta t=0.001\) |
| Pressure-law exponent | \(m=80\) |
| gPC polynomial degree | \(K=4\) |
| Reference method | 16-node Gauss–Legendre stochastic collocation |
| Initial core and rim densities | \(0.15-0.1z\) and \(0.35+0.1z\) |
| Initial radial parameters | \((0.4,0.5,0.6,0.7)(1+0.3z)\) |
| Growth coefficient | \(G_0=0.5(1-0.1z)\) |

The piecewise initial density is defined in equation (29). The experiment examines density means, standard deviations, and spatial slices. The paper also reports errors with increasing polynomial degree; this repository focuses on the selected concave experiment rather than the full parameter study.

## Running the Code

In MATLAB, set the current folder to the repository root.

```matlab
run_experiment                     % Recompute SG and all SC trajectories at three times
run_experiment('cached_reference')  % Recompute SG using the saved SC reference
run_experiment('reference_short')   % Recompute SG and SC at the first output time
plot_saved                         % Plot the supplied reference datasets
run_smoke                          % Run a short solver and plotting check
```

`run_experiment` performs initial projection, time integration, statistical reconstruction, and plotting. Inputs are loaded from `data/concave_compute.mat`. Generated results are written to `results/` without overwriting the supplied datasets.

`plot_saved` generates figures directly from the saved experimental data. `run_smoke` runs SG and all 16 SC nodes to \(T=0.01\), then checks the saved-data plotting workflow.

## Repository Contents

| Path | Purpose |
|---|---|
| `run_experiment.m` | Main computation and statistical post-processing |
| `run_smoke.m` | Short execution check |
| `plot_saved.m` | Plotting from saved experimental data |
| `solver/solve_gpc_sg.m` | Coupled gPC-SG solver |
| `solver/solve_sc.m` | Deterministic solver for each collocation node |
| `solver/gpc_growth.m` | Growth-term projection for SG |
| `solver/sc_growth.m` | Growth-term evaluation for SC |
| `solver/lgwt.m` | Gauss–Legendre nodes and weights |
| `data/concave_compute.mat` | Computational inputs and reference statistics |
| `data/gpc_concave_paper.mat` | Saved data for the selected paper figures |
| `TESTS.md` | Execution checks and test details |

The main entry point calls the SG and SC solvers and their growth-term routines. Auxiliary flux and matrix operations remain within the corresponding solver files.

The experiment originates from the concave branch of `test1111.m`; the saved plotting dataset originates from `1127concave_max_vivo.mat`.

Keep `solver/` and `data/` in their existing locations. Entry points resolve paths relative to the repository, and the generated `results/` directory is excluded by `.gitignore`.
