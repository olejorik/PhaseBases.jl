# PhaseBases — TODO

## Generalize coordinate transforms for Zernike basis

**Current state:** `ZernikeBW(dom, d, order; coordmap=identity)` accepts an arbitrary
function `coordmap((x,y)) -> (u,v)` that transforms coordinates before evaluating
Zernike polynomials. Pre-defined maps: `rot90ccw`, `rot90cw`, `rot180`, `flipx`, `flipy`.

**Limitation:** `coordmap` is a plain function — no inverse, no composition, no
reflection in the `ZernikeBW` type itself. There is no way to ask "what transform was
used to build this basis?".

**Desired improvement:**
- Introduce a proper `CoordMap` type (or reuse `AffineDomain2D` from SampledDomains
  once it exists) so the transform is inspectable and composable.
- Consider storing the `coordmap` inside the `ZernikeBW` struct so it is always
  self-describing.
- Provide a general `AffineCoordMap(A::SMatrix{2,2})` for arbitrary linear transforms
  (rotation by arbitrary angle, anisotropic scaling, shear).

**Related:** `AffineDomain2D` in SampledDomains (see that package's TODO).

## Parametrize `Basis`/`ZernikeBW` struct fields to eliminate allocations in `decompose!`

**Current state:** `decompose!`/`decompose`/`allinners!` for the general
`AbstractBasis` fallback now dispatch via a `BasisLayoutStyle` trait
(`Indexed`/`Masked`/`Unsupported`, see `basislayoutstyle` in `types.jl`) instead of
runtime `hasproperty` checks. This is cleaner but did **not** eliminate the
allocation in `decompose!`/`project!`/`residual!` for `Basis` (~624 bytes/call
measured).

**Root cause:** `fieldtypes(Basis)` are not concrete:
`(VectorOfArray, VectorOfArray, Array, Vector{<:CartesianIndex}, Vector)` — no type
parameters on `VectorOfArray`/`Array`/`Vector`, and `Vector{<:CartesianIndex}` is a
`UnionAll`, not a concrete type. Any access to `b.dualelements`/`b.indexes` is
therefore inferred as abstract, forcing the compiler to box every scalar produced
inside the `decompose!` loop (`inner_indexed` result) — independent of the trait
dispatch, since the trait only picks which method runs, not the field types.
`ZernikeBW` likely has the same issue (same field pattern) and should be checked too.

**Suggested fix:** parametrize the struct with concrete field types, e.g.:

```julia
struct Basis{TE<:VectorOfArray,TD<:VectorOfArray,TA<:AbstractArray,
             TI<:AbstractVector{<:CartesianIndex},TN<:AbstractVector} <: AbstractBasis
    elements::TE
    dualelements::TD
    ap::TA
    indexes::TI
    norms::TN
end
```

Then re-verify `decompose!`/`project!`/`residual!` allocate 0 bytes via `@allocated`
inside a local function (not at global/top-level scope, which itself causes
unrelated spurious allocations). Apply the same treatment to `ZernikeBW` if it
shows the same issue.

