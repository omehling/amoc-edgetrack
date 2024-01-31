This repository provides Julia source code for the paper

O. Mehling, R. Börner, V. Lucarini (2024): **Limits to predictability of the asymptotic state of the Atlantic Meridional Overturning Circulation in a conceptual climate model.** *Physica D: Nonlinear Phenomena* 459, 134043, https://doi.org/10.1016/j.physd.2023.134043

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10370900.svg)](https://doi.org/10.5281/zenodo.10370900)

## Overview

The coupled conceptual Stommel–L84 model as proposed by Gottwald ([2021](https://doi.org/10.1007/s00382-020-05476-z)) is implemented in `model_definitions.jl`.

The other files provide methods and scripts to reproduce different results from the paper:
* `edgetrack_L84stommel_sigma.jl` runs the edge tracking algorithm and produces a pseudo-trajectory of the chaotic saddle for a given value of $\sigma_0$
* `escape_rates.jl` determines the escape rate from the saddle given a pseudo-trajectory of the saddle
* `lyapunov_saddle.jl` determines the Lyapunov spectrum given a pseudo-trajectory of the saddle
* `transient_lifetime_*.jl` calculate the lifetimes of long transients outside of the bistable regime
* `basin_boundary_section.jl` samples different phase space sections through the fractal basin boundary to compute its box-counting dimension

## Installation

Code was written for and run with **Julia 1.8.1** using [`DynamicalSystems.jl`](https://github.com/JuliaDynamics/DynamicalSystems.jl) version 2.3. It may not work with other versions, and specifically it is not compatible with DynamicalSystems 3.X.

For reproducibility, we provide the files `Project.toml` and `Manifest.toml` to replicate the Julia environment via `Pkg.instantiate`, and recommend [`juliaup`](https://github.com/JuliaLang/juliaup) to set the Julia version of the working directory:

```
$ git clone https://github.com/omehling/amoc-edgetrack.git
$ cd amoc-edgetrack

$ juliaup add 1.8.1            # Install Julia version 1.8.1 if not available already
$ juliaup override set 1.8.1   # Set Julia version for this directory

$ julia
julia> ]
pkg> activate .
pkg> instantiate
```

## Notes

The edge tracking algorithm given in `edgetrack.jl` is superseded by the [`edgetracking`](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/attractors/stable/basins/#Attractors.edgetracking) method of [`Attractors.jl`](https://github.com/JuliaDynamics/Attractors.jl).

It is **strongly recommended** that new code uses this new implementation instead of the one provided here.
