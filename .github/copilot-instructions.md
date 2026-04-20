# PhaseBases.jl — Copilot Context

## What this package is

`PhaseBases` provides basis representations for optical phase aberrations, used in wavefront sensing, phase retrieval, and adaptive optics. It is a foundational dependency of `PhaseRetrieval` and `Feedback14AMI`.

## Key types and functions

| Symbol | Purpose |
|---|---|
| `Basis` / `AbstractBasis` | Interface: collection of 2D arrays with aperture weighting |
| `ZernikeBW` | Zernike basis, Born & Wolf indexing, normalized on unit disk |
| `ZernikeBWSparse` | Sparse variant for large grids |
| `PixelBasis` | Zonal/pixel-by-pixel basis |
| `ShiftedBasis` | Basis shifted to a sub-aperture |
| `ModalPhase` | Phase as linear combination of basis elements; holds `coef` + `basis` |
| `ZonalPhase` | Pixel-map phase representation |
| `SymbolicZernikePhase` | Coefficient-only description; deferred grid discretization |
| `compose` / `compose!` | Reconstruct phase array from coefficients |
| `decompose` / `decompose!` | Project phase array onto basis (inner products) |
| `allinners!` | Batch inner products; used in gradient projection |
| `dualbasis` | Compute dual basis (inverse Gram matrix) |

## Zernike ordering conventions

| Convention | Tag | First index |
|---|---|---|
| OSA / ANSI | `OSA` | 0 |
| Noll | `Noll` | 1 |
| Fringe (U. Arizona) | `Fringe` | 1 |
| Born & Wolf | `BornWolf` | 1 |

## Array axis convention

All 2D arrays follow **row = y, column = x**: `A[j, i] = f(x[i], y[j])`. For Makie `heatmap`, arrays must be transposed for correct spatial orientation.

## Relationships

- Used by: `PhaseRetrieval`, `Feedback14AMI` (gradient projection onto Zernike basis)
- Depends on: `SampledDomains` (grid description), `RecursiveArrayTools`

## Source layout

```
src/
    PhaseBases.jl   ← module entry, exports
    types.jl        ← Phase, ModalPhase, ZonalPhase, abstract types
    Basis.jl        ← Basis interface, compose/decompose, allinners!
    Zernike/        ← ZernikeBW, symbolic Zernike, ordering conversions
examples/           ← standalone runnable examples
```
