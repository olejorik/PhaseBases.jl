# About PhaseBases.jl

`PhaseBases` provides basis representations for optical phase aberrations,
intended for use in wavefront sensing, phase retrieval, and adaptive optics
simulations.

## Key design ideas

- **Bases are vector-like**: a basis is a collection of arrays; indexing,
  composition, and decomposition all follow a uniform interface.
- **Phase types are separate from arrays**: `ModalPhase` and `ZonalPhase`
  hold a lightweight description of a wavefront; `collect` materializes them
  to an array when needed.
- **Symbolic coefficients first**: `SymbolicZernikePhase` lets you specify a
  wavefront purely in terms of Zernike coefficients and an ordering convention,
  deferring grid discretization to a later step.

## Supported Zernike conventions

| Convention | Type tag | First index |
|---|---|---|
| OSA / ANSI | `OSA` | 0 |
| Noll | `Noll` | 1 |
| Fringe (University of Arizona) | `Fringe` | 1 |
| Mizer | `Mizer` | 1 |

## Array axis convention

Currently, all 2D arrays produced by this package (Zernike mode images, apertures, phase maps) follow the **row = y, column = x** convention:

```
A[j, i] = f(x[i], y[j])
```

- `i` indexes **x** (`xrange[i]`), the horizontal / column direction — second array index.
- `j` indexes **y** (`yrange[j]`), the vertical / row direction — first array index.

This matches how Julia (and `Images.jl`) loads and displays images: a sensor image that is `W` pixels wide and `H` pixels tall is a `H × W` matrix, with `img[row, col]` addressing pixel at vertical position `row` and horizontal position `col`. The REPL display is therefore spatially correct — what you see printed is what you get in the optical plane.


To display the calculated Zernike polynomials with Makie's `heatmap` function in correct orientation, the arrays should be transposed

## Related packages

- [`PhaseUtils`](https://github.com/olejorik/PhaseUtils.jl) — phase unwrapping and windowing utilities
- [`PhaseRetrieval`](https://github.com/olejorik/PhaseRetrieval.jl) — phase retrieval algorithms
- [`PhasePlots`](https://github.com/olejorik/PhasePlots.jl) — visualization helpers
