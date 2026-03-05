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

## Related packages

- [`PhaseUtils`](https://github.com/olejorik/PhaseUtils.jl) — phase unwrapping and windowing utilities
- [`PhaseRetrieval`](https://github.com/olejorik/PhaseRetrieval.jl) — phase retrieval algorithms
- [`PhasePlots`](https://github.com/olejorik/PhasePlots.jl) — visualization helpers
