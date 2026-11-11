# porous-medium-tumor-gpc-sg-uq
This repository provides a selected numerical experiment accompanying:

**Ning Jiang, Liu Liu, and Huimin Yu. “Efficient stochastic asymptotic-preserving scheme for tumor growth models with uncertain parameters.” Journal of Computational and Applied Mathematics, 473 (2026), 116881.**  
https://doi.org/10.1016/j.cam.2025.116881

The paper studies uncertainty quantification for tumor growth models governed by porous-medium-type equations coupled with nutrient-dependent growth. Uncertainties in the initial cell density, tumor support, and growth parameters propagate through the evolving tumor interface and influence the resulting density statistics. The proposed method combines generalized polynomial chaos with stochastic Galerkin projection in a stochastic asymptotic-preserving framework. The analysis connects the stochastic Galerkin system to the corresponding Hele–Shaw limit as the pressure-law exponent increases.

This repository focuses on the **concave initial-density experiment in Test I(a)**. The experiment considers a two-dimensional spatial domain and one uniformly distributed random variable affecting the initial density, radial support, and growth coefficient. It uses a fourth-order gPC approximation with pressure-law exponent \(m=80\), and compares the computed mean and standard deviation with a 16-node Gauss–Legendre stochastic collocation reference.

The implementation includes initial-data projection, the coupled gPC-SG solver, the stochastic collocation solver, statistical post-processing, and plotting routines. Both recomputation from the supplied numerical inputs and direct plotting of the included reference datasets are supported. The selected experiment illustrates how random inputs affect tumor-density statistics and assesses the accuracy of the gPC-SG approximation against stochastic collocation.
