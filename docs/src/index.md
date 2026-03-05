```@meta
CurrentModule = PhaseBases
```

# PhaseBases.jl

Basis representations for optical phase aberrations.

`PhaseBases` provides types and functions for working with sets of functions
(bases) defined on a computational aperture, with a focus on Zernike polynomials
and modal phase representations.

## Highlights

- [`ZernikeBW`](@ref) — disk-based Zernike polynomial basis (Born & Wolf normalization)
- [`SymbolicZernikePhase`](@ref) — coefficients-only phase, independent of any computational grid
- [`ModalPhase`](@ref) — phase as a linear combination of basis functions
- [`ZonalPhase`](@ref) — phase as a pixel-value array
- Index conversion utilities: `fringe_j_to_nm`, `noll_j_to_nm`, `nm_to_osa_j`, …

## Quick start

```julia
using PhaseBases

# Create a Zernike basis on a 128×128 grid, radial orders 0..6
zbas = ZernikeBW(128, 6)

# Describe a wavefront symbolically (Fringe coefficients from an interferometer)
wf = SymbolicZernikePhase([0.02, -0.1, 0.05, 0.8, -0.3, 0.0, 0.1, 0.0, 0.5], Fringe)

# Evaluate on the grid
arr = collect(wf, zbas)
```

## Index

```@index
```
