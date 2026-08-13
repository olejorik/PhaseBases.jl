# Changelog

All notable changes to PhaseBases.jl are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
- `zernike_basis(indices, ordering, dom, d; ...)` — build a `Basis` restricted to a specific set of Zernike modes (by single index, in any `ZernikeOrdering`), without paying the `pinv` cost of a full triangular `ZernikeBW` basis
- `zernike_basis(nterms, ordering, dom, d; ...)` — dense convenience for the first `nterms` consecutive indices in a given ordering (e.g. the first 36 Fringe terms)
- `zernike_basis(ph::SymbolicZernikePhase, dom, d; ...)` — build the basis matching a symbolic phase's own modes exactly
- `ModalPhase(ph::SymbolicZernikePhase, dom, d; ...)` — convert a symbolic phase directly to a `ModalPhase` with no zero-padding, using the modes-matching basis above
- `CHANGELOG.md` — started tracking changes following Keep a Changelog format
- `.github/copilot-instructions.md` — added project context for GitHub Copilot (key types, exported functions, relationships, source layout)

