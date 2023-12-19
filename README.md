This repository provides Julia source code for the paper

O. Mehling, R. Börner, V. Lucarini (2023): **Limits to predictability of the asymptotic state of the Atlantic Meridional Overturning Circulation in a conceptual climate model** [arXiv:2308.16251](http://arxiv.org/abs/2308.16251)

---

The coupled conceptual Stommel–L84 model as proposed by Gottwald ([2021](https://doi.org/10.1007/s00382-020-05476-z)) is implemented in `model_definitions.jl`.

The other files provide methods and scripts to reproduce different results from the paper:
* `edgetrack_L84stommel_sigma.jl` runs the edge tracking algorithm and produces a pseudo-trajectory of the chaotic saddle for a given value of $\sigma_0$
* `escape_rates.jl` determines the escape rate from the saddle given a pseudo-trajectory of the saddle
* `lyapunov_saddle.jl` determines the Lyapunov spectrum given a pseudo-trajectory of the saddle
* `transient_lifetime_*.jl` calculate the lifetimes of long transients outside of the bistable regime
* `basin_boundary_section.jl` samples different phase space sections through the fractal basin boundary to compute its box-counting dimension

All scripts build on the [`DynamicalSystems.jl`](https://github.com/JuliaDynamics/DynamicalSystems.jl) library **version 2.3**.

**NB:** The edge tracking algorithm given in `edgetrack.jl` is superseded by a new implementation of the `edgetracking` method to be released as part of DynamicalSystems.jl: https://github.com/reykboerner/Attractors.jl/blob/85fb37068bd31cc97d4091e6452edb7c4c255af9/src/boundaries/edgetracking.jl <br />
It is strongly recommended that new code uses this implementation instead of the one provided here.
