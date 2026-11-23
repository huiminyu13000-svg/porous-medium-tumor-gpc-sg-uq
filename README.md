# Porous-Medium Tumor Growth UQ with gPC Stochastic Galerkin

This repository accompanies N. Jiang, L. Liu, and H. Yu, *Efficient stochastic asymptotic-preserving scheme for tumor growth models with uncertain parameters*, Journal of Computational and Applied Mathematics **473** (2026), 116881. It contains one complete experiment: the single-input concave initial density in **Test I(a)**, corresponding to the left-hand example on slide 7 of the defense presentation.

The program generates the initial condition, projects it onto a Legendre basis, computes the stochastic Galerkin (SG) and stochastic collocation (SC) solutions, and displays density means and standard deviations. All inputs are generated from formulas; no input MAT files are required.

## Run

Set the MATLAB current folder to this repository root:

```matlab
run_experiment          % SG and all 16 SC nodes at T = 0.01, 1, 2
run_experiment('short') % The same solvers at T = 0.01
```

Figures are displayed and saved as PNG files in `results/`, together with the calculated statistics, inputs, settings, and timings. The plotting functions are at the end of `run_experiment.m`. MATLAB or GNU Octave can be used; no additional MATLAB toolbox is needed. Octave needs an available graphics toolkit for plotting.

## Experiment settings

| Quantity | Setting |
|---|---|
| Random input | One variable, $z\sim\mathcal U[-1,1]$ |
| Spatial domain and spacing | $[-2.2,2.2]^2$, $\Delta x=\Delta y=0.1$ |
| Time step and pressure exponent | $\Delta t=0.001$, $m=80$ |
| SG expansion | Four orthonormal Legendre modes, degrees 0 through 3 |
| SC reference | 16-point Gauss–Legendre quadrature, probability weights summing to 1 |
| Growth projection | 32-point Gauss–Legendre quadrature |
| Output times | 0.01, 1, 2; spatial fields at $T=1$ |
| Statistical slices | $y=-0.3$ |

With $r=\sqrt{x^2+y^2}$, define

$$
\rho_c=0.15-0.1z,\qquad \rho_e=0.35+0.1z,\qquad
(r_1,r_2,r_3,r_4)=(0.4,0.5,0.6,0.7)(1+0.3z).
$$

The initial density is $\rho_c$ for $r\le r_1$, interpolates linearly from $\rho_c$ to $\rho_e$ on $(r_1,r_2]$, equals $\rho_e$ on $(r_2,r_3]$, and decreases linearly to zero on $(r_3,r_4]$. It is zero for $r>r_4$.

The growth field is $h=G_0c$, with $G_0=0.5(1-0.1z)$. The radial in-vivo nutrient field uses interior consumption coefficient 1 and exterior supply $c_B=1$. Its radius is obtained from the current positive support (the mean density for SG, each individual solution for SC). The SG solver retains the supplied coefficientwise pressure discretization and coupled growth projection. SG standard deviation is the square root of the sum of squared nonconstant coefficients; SC statistics use the quadrature weights.

## Code structure

| File | Purpose |
|---|---|
| `experiment_config.m` | Experiment settings |
| `run_experiment.m` | Input generation, projection, SG/SC calculation, statistics, and plots |
| `solver/solve_gpc_sg.m`, `solver/solve_sc.m` | Front-capturing solvers, including flux and matrix operations |
| `solver/initial_density.m`, `solver/nutrient_field.m` | Initial condition and nutrient model |
| `solver/gpc_growth.m`, `solver/sc_growth.m` | Growth terms |
| `solver/legendre_basis.m`, `solver/lgwt.m` | Orthonormal basis and quadrature |

The final time step is shortened if needed to reach the requested output time exactly. This repository contains the selected single-input experiment; the two-input and multi-species experiments are outside its scope. Generated `results/` files are excluded from Git.
