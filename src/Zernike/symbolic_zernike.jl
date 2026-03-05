# ===============================================================================
# Symbolic Zernike Phase — basis-free modal representation
# ===============================================================================
#
# This file defines `SymbolicZernikePhase`, a phase type that holds Zernike
# coefficients in a specific ordering convention (Fringe, Noll, OSA, Mizer)
# without reference to any discretized grid or concrete basis object.
#
# Conversion to a concrete array or ModalPhase happens separately, when a
# grid (CartesianDomain2D + diameter) or an existing ZernikeBW basis is provided.

export SymbolicZernikePhase, ZernikeOrdering, ZernikeNormalization
export FringeOrdering, NollOrdering, OSAOrdering, MizerOrdering
export reorder, to_nm


# -------------------------------------------------------------------------------
# Ordering and normalization type tags (open for extension via subtyping)
# -------------------------------------------------------------------------------

"""
    ZernikeOrdering

Abstract type for Zernike polynomial single-index ordering conventions.
Concrete subtypes (singleton structs) select a convention:
- `FringeOrdering` / `Fringe` — University of Arizona Fringe convention (1-based)
- `NollOrdering`   / `Noll`  — Noll's convention (1-based)
- `OSAOrdering`    / `OSA`   — OSA/ANSI standard (0-based)
- `MizerOrdering`  / `Mizer` — Mizer internal convention (1-based)

New conventions can be added by defining a new subtype and the four methods:
`j_to_nm(ordering, j)`, `nm_to_j(ordering, nm)`, and `j_first(ordering)`.
"""
abstract type ZernikeOrdering end

struct FringeOrdering <: ZernikeOrdering end
struct NollOrdering   <: ZernikeOrdering end
struct OSAOrdering    <: ZernikeOrdering end
struct MizerOrdering  <: ZernikeOrdering end

"""
    ZernikeNormalization

Abstract type for Zernike polynomial normalization conventions.
- `BWNormalization`  / `BornWolf` — unit amplitude at aperture edge (classical)
- `RMSNormalization` / `RMSNorm`  — unit RMS over unit disk (OSA standard)
"""
abstract type ZernikeNormalization end

struct BWNormalization  <: ZernikeNormalization end
struct RMSNormalization <: ZernikeNormalization end

## Convenience singleton constants — use these in practice
const Fringe   = FringeOrdering()
const Noll     = NollOrdering()
const OSA      = OSAOrdering()
const Mizer    = MizerOrdering()
const BornWolf = BWNormalization()
const RMSNorm  = RMSNormalization()

export Fringe, Noll, OSA, Mizer, BornWolf, RMSNorm


# -------------------------------------------------------------------------------
# Ordering dispatch: single index ↔ (n, m)
# -------------------------------------------------------------------------------

"""
    j_to_nm(ordering::ZernikeOrdering, j::Int) -> NamedTuple{(:n, :m), ...}

Convert a single Zernike index `j` in the given `ordering` to the double (n, m) index.
"""
j_to_nm(::FringeOrdering, j::Int) = fringe_j_to_nm(j)
j_to_nm(::NollOrdering,   j::Int) = noll_j_to_nm(j)
j_to_nm(::OSAOrdering,    j::Int) = osa_j_to_nm(j)
j_to_nm(::MizerOrdering,  j::Int) = mizer_j_to_nm(j)

"""
    nm_to_j(ordering::ZernikeOrdering, nm::NamedTuple) -> Int

Convert double (n, m) index to a single index in the given `ordering`.
"""
nm_to_j(::FringeOrdering, nm) = nm_to_fringe_j(nm)
nm_to_j(::NollOrdering,   nm) = nm_to_noll_j(nm)
nm_to_j(::OSAOrdering,    nm) = nm_to_osa_j(nm)
nm_to_j(::MizerOrdering,  nm) = nm_to_mizer_j(nm)

"""
    j_first(ordering::ZernikeOrdering) -> Int

Return the first (lowest) single index value for the given ordering.
"""
j_first(::OSAOrdering)   = 0  ## 0-based
j_first(::FringeOrdering) = 1
j_first(::NollOrdering)   = 1
j_first(::MizerOrdering)  = 1


# -------------------------------------------------------------------------------
# SymbolicZernikePhase
# -------------------------------------------------------------------------------

"""
    SymbolicZernikePhase{T<:Real, TO<:ZernikeOrdering, TN<:ZernikeNormalization}

A phase aberration described by Zernike polynomial coefficients without
reference to any discretized grid or concrete basis object.

Stores coefficients sparsely as (index → value) pairs in a specific ordering
convention. Can be converted to a concrete `ModalPhase` or phase array once a
computational grid is provided.

# Constructors

    SymbolicZernikePhase(indices, coef, ordering [, normalization])

Sparse: coefficients for the Zernike polynomials at the given `indices`.

    SymbolicZernikePhase(coef, ordering [, normalization])

Dense: coefficients for *consecutive* indices starting from `j_first(ordering)`.

# Examples
```julia
## Single defocus term, Noll index 4
σ = SymbolicZernikePhase([4], [1.0], Noll)

## First 15 Fringe coefficients (dense, from index 1)
σ = SymbolicZernikePhase([0.1, 0.0, -0.2, 0.5, 0.0, 0.0, 0.0, 0.0, 0.8], Fringe)

## Specific Fringe terms only (sparse)
σ = SymbolicZernikePhase([1, 4, 9], [0.0, 1.5, 0.3], Fringe)
```
"""
struct SymbolicZernikePhase{T<:Real, TO<:ZernikeOrdering, TN<:ZernikeNormalization} <: AbstractModalPhase
    indices::Vector{Int}
    coef::Vector{T}
    ordering::TO
    normalization::TN

    function SymbolicZernikePhase(
        indices::AbstractVector{<:Integer},
        coef::AbstractVector{T},
        ordering::TO,
        normalization::TN=BornWolf,
    ) where {T<:Real, TO<:ZernikeOrdering, TN<:ZernikeNormalization}
        length(indices) == length(coef) ||
            throw(ArgumentError("indices and coef must have the same length (got $(length(indices)) and $(length(coef)))"))
        allunique(indices) ||
            throw(ArgumentError("Zernike indices must be unique"))
        return new{T,TO,TN}(collect(Int, indices), collect(T, coef), ordering, normalization)
    end
end

## Dense constructor: coefficients for consecutive indices starting from j_first(ordering)
function SymbolicZernikePhase(
    coef::AbstractVector{T},
    ordering::TO,
    normalization::TN=BornWolf,
) where {T<:Real, TO<:ZernikeOrdering, TN<:ZernikeNormalization}
    j0 = j_first(ordering)
    indices = collect(j0:(j0 + length(coef) - 1))
    return SymbolicZernikePhase(indices, coef, ordering, normalization)
end


# -------------------------------------------------------------------------------
# Phase interface
# -------------------------------------------------------------------------------

coefficients(ph::SymbolicZernikePhase) = ph.coef

function coefficients!(ph::SymbolicZernikePhase, coef)
    length(coef) == length(ph.coef) ||
        throw(ArgumentError("New coefficient vector must have the same length as the existing one"))
    ph.coef .= coef
    return ph
end

copy(ph::SymbolicZernikePhase) =
    SymbolicZernikePhase(copy(ph.indices), copy(ph.coef), ph.ordering, ph.normalization)

function -(ph::SymbolicZernikePhase)
    ph1 = copy(ph)
    ph1.coef .*= -1
    return ph1
end

function *(c::Real, ph::SymbolicZernikePhase)
    ph1 = copy(ph)
    ph1.coef .*= c
    return ph1
end

*(ph::SymbolicZernikePhase, c::Real) = c * ph

function show(io::IO, ph::SymbolicZernikePhase{T,TO,TN}) where {T,TO,TN}
    nz = count(!iszero, ph.coef)
    n  = length(ph.coef)
    print(io, "SymbolicZernikePhase{$T} — $nz/$n nonzero term(s), ordering: $(nameof(TO)), normalization: $(nameof(TN))")
end


# -------------------------------------------------------------------------------
# Index conversion
# -------------------------------------------------------------------------------

"""
    to_nm(ph::SymbolicZernikePhase) -> Vector{NamedTuple{(:n,:m)}}

Return all coefficient indices converted to (n, m) double-index form.
"""
to_nm(ph::SymbolicZernikePhase) = [j_to_nm(ph.ordering, j) for j in ph.indices]


"""
    reorder(ph::SymbolicZernikePhase, new_ordering::ZernikeOrdering) -> SymbolicZernikePhase

Convert a symbolic Zernike phase to a different single-index ordering convention.
Coefficients are unchanged; only the indices are remapped via (n, m).
"""
function reorder(ph::SymbolicZernikePhase, new_ordering::ZernikeOrdering)
    typeof(ph.ordering) == typeof(new_ordering) && return copy(ph)
    nm_pairs    = to_nm(ph)
    new_indices = [nm_to_j(new_ordering, nm) for nm in nm_pairs]
    return SymbolicZernikePhase(new_indices, copy(ph.coef), new_ordering, ph.normalization)
end


# -------------------------------------------------------------------------------
# Arithmetic
# -------------------------------------------------------------------------------

"""
    ph1 + ph2 :: SymbolicZernikePhase

Add two symbolic Zernike phases. If they use different orderings, `ph2` is
reindexed to `ph1`'s convention before addition. Adding phases with different
normalizations is an error — convert explicitly first.
"""
function +(ph1::SymbolicZernikePhase, ph2::SymbolicZernikePhase)
    typeof(ph1.normalization) == typeof(ph2.normalization) ||
        error("Cannot add SymbolicZernikePhases with different normalizations " *
              "($(nameof(typeof(ph1.normalization))) vs $(nameof(typeof(ph2.normalization)))). Convert first.")

    ## Bring ph2 into ph1's ordering
    ph2c = reorder(ph2, ph1.ordering)

    ## Merge index sets (preserve ph1 ordering, append new ones from ph2)
    idx1 = Dict(j => i for (i, j) in enumerate(ph1.indices))
    T    = promote_type(eltype(ph1.coef), eltype(ph2c.coef))

    ## Start from ph1's full coefficient vector
    merged_indices = copy(ph1.indices)
    merged_coef    = Vector{T}(ph1.coef)

    for (k, j) in enumerate(ph2c.indices)
        if haskey(idx1, j)
            merged_coef[idx1[j]] += ph2c.coef[k]
        else
            push!(merged_indices, j)
            push!(merged_coef, T(ph2c.coef[k]))
        end
    end

    return SymbolicZernikePhase(merged_indices, merged_coef, ph1.ordering, ph1.normalization)
end

-(ph1::SymbolicZernikePhase, ph2::SymbolicZernikePhase) = ph1 + (-ph2)


# -------------------------------------------------------------------------------
# Materialization
# -------------------------------------------------------------------------------

"""
    collect(ph::SymbolicZernikePhase, basis::ZernikeBW) -> Array

Evaluate the symbolic phase on the grid defined by `basis`. Indices are
converted via (n, m) to the OSA ordering used internally by `ZernikeBW`
(1-based array position = OSA index + 1).

Note: normalization conversion (BornWolf ↔ RMSNorm) is not yet implemented.
"""
function collect(ph::SymbolicZernikePhase, basis::ZernikeBW)
    nm_pairs    = to_nm(ph)
    osa_indices = [nm_to_osa_j(nm) + 1 for nm in nm_pairs]   ## 1-based array position

    full_coef = zeros(eltype(ph.coef), length(basis))
    for (k, idx) in enumerate(osa_indices)
        1 ≤ idx ≤ length(basis) ||
            throw(ArgumentError(
                "Zernike term $(nm_pairs[k]) (array index $idx) " *
                "exceeds basis length $(length(basis))"
            ))
        full_coef[idx] += ph.coef[k]
    end

    return compose(basis, full_coef)
end

"""
    collect(ph::SymbolicZernikePhase, dom::CartesianDomain2D, diameter::Real) -> Array

Convenience method: construct a `ZernikeBW` basis of sufficient order from
the domain and circular aperture diameter, then evaluate the phase.
"""
function collect(ph::SymbolicZernikePhase, dom::CartesianDomain2D, diameter::Real)
    nm_pairs = to_nm(ph)
    maxord   = maximum(nm.n for nm in nm_pairs)
    basis    = ZernikeBW(dom, diameter, maxord)
    return collect(ph, basis)
end

"""
    ModalPhase(ph::SymbolicZernikePhase, basis::ZernikeBW) -> ModalPhase

Convert a symbolic Zernike phase to a concrete `ModalPhase` in `basis`.
"""
function ModalPhase(ph::SymbolicZernikePhase, basis::ZernikeBW)
    nm_pairs    = to_nm(ph)
    osa_indices = [nm_to_osa_j(nm) + 1 for nm in nm_pairs]

    full_coef = zeros(eltype(ph.coef), length(basis))
    for (k, idx) in enumerate(osa_indices)
        1 ≤ idx ≤ length(basis) ||
            throw(ArgumentError(
                "Zernike term $(nm_pairs[k]) (array index $idx) " *
                "exceeds basis length $(length(basis))"
            ))
        full_coef[idx] += ph.coef[k]
    end

    return ModalPhase(full_coef, basis)
end
