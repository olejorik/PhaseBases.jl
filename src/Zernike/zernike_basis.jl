# ===============================================================================
# zernike_basis — construct a `Basis` restricted to a specific set of Zernike modes
# ===============================================================================
#
# Unlike `ZernikeBW(dom, d, maxorder)`, which builds the full triangular set of
# modes up to `maxorder` and computes the dual (pseudo-inverse) basis over ALL of
# them, these constructors select only the requested modes *before* computing the
# dual basis, so the cost of the `pinv` scales with the number of modes actually
# requested, not with the size of the triangular basis needed to reach the same
# maximum radial order.

export zernike_basis

"""
    zernike_basis(indices, ordering::ZernikeOrdering, dom::CartesianDomain2D, d::Real;
                  normalization::ZernikeNormalization=BornWolf,
                  fftshifted::Bool=false, coordmap=identity) -> Basis

Build a `Basis` containing only the Zernike polynomials at `indices`, given in the
single-index convention `ordering` (e.g. [`Fringe`](@ref), [`Noll`](@ref),
[`OSA`](@ref), [`Mizer`](@ref)), sampled on `dom` and restricted to a circular
aperture of diameter `d`.

Elements in the returned `Basis` follow the order of `indices`; `indices` must be
unique (duplicates raise an `ArgumentError`).

The dual (pseudo-inverse) basis is computed only over the requested modes, so this
is substantially cheaper than building a full `ZernikeBW(dom, d, maxorder)` and
discarding unused elements — the `pinv` cost scales with `length(indices)`, not
with the size of the full triangular basis needed to reach the same maximum radial
order. Note that generating the raw Zernike values themselves is not reduced (the
recurrence used by [`zernike`](@ref) computes all lower orders together), only the
subsequent dual-basis computation.

Only `BornWolf` normalization is currently supported; passing `RMSNorm` raises an
error.

Set `fftshifted=true` to apply `ifftshift` to the elements and aperture, matching
the FFT-centered pupil convention used e.g. by
`ZernikeBW(dom, d, maxorder; fftshifted=true)`.

Note: unlike `ZernikeBW`, the returned `Basis` has no `mask` field for NaN-masked
plotting; use `aperture(basis)` for a 0/1 mask instead.

# Examples
```julia
julia> dom = PhaseBases.CartesianDomain2D(-1:0.05:1, -1:0.05:1);

julia> bas = zernike_basis(1:9, Fringe, dom, 1.6);   # first 9 Fringe terms only

julia> length(bas)
9
```

See also [`ZernikeBW`](@ref), [`SymbolicZernikePhase`](@ref).
"""
function zernike_basis(
    indices::AbstractVector{<:Integer},
    ordering::ZernikeOrdering,
    dom::CartesianDomain2D,
    d::Real;
    normalization::ZernikeNormalization=BornWolf,
    fftshifted::Bool=false,
    coordmap=identity,
)
    normalization === BornWolf ||
        error("zernike_basis currently supports only BornWolf normalization")
    allunique(indices) || throw(ArgumentError("indices must be unique"))

    nm_pairs = [j_to_nm(ordering, j) for j in indices]
    maxorder = maximum(nm.n for nm in nm_pairs)

    ztable = makezerniketable(dom, maxorder, d / 2; coordmap)
    osa_positions = [nm_to_osa_j(nm) + 1 for nm in nm_pairs]   ## 1-based array position

    els = [ztable[pos] for pos in osa_positions]
    ap, _ = aperture(dom, d)

    if fftshifted
        els = [ifftshift(e) for e in els]
        ap = ifftshift(ap)
    end

    idx = findall(!iszero, ap)
    return Basis(els, idx)
end

"""
    zernike_basis(nterms::Integer, ordering::ZernikeOrdering,
                  dom::CartesianDomain2D, d::Real; kwargs...) -> Basis

Dense convenience: build a `Basis` from the first `nterms` consecutive indices in
`ordering`, starting at `j_first(ordering)` — matching the dense constructor of
[`SymbolicZernikePhase`](@ref).

For `Fringe`, the natural truncation points are the ring boundaries
`nterms = p^2` (9, 16, 25, 36, ...); other values of `nterms` are still valid,
just not a complete Fringe ring.
"""
function zernike_basis(
    nterms::Integer, ordering::ZernikeOrdering, dom::CartesianDomain2D, d::Real; kwargs...
)
    j0 = j_first(ordering)
    return zernike_basis(j0:(j0 + nterms - 1), ordering, dom, d; kwargs...)
end

"""
    zernike_basis(ph::SymbolicZernikePhase, dom::CartesianDomain2D, d::Real;
                  fftshifted::Bool=false, coordmap=identity) -> Basis

Build a `Basis` containing exactly the modes present in `ph` (its own `indices`, in
its own `ordering` and `normalization`), suitable for constructing a `ModalPhase`
matching `ph` with no zero-padding.
"""
function zernike_basis(
    ph::SymbolicZernikePhase,
    dom::CartesianDomain2D,
    d::Real;
    fftshifted::Bool=false,
    coordmap=identity,
)
    return zernike_basis(
        ph.indices,
        ph.ordering,
        dom,
        d;
        normalization=ph.normalization,
        fftshifted=fftshifted,
        coordmap=coordmap,
    )
end

"""
    ModalPhase(ph::SymbolicZernikePhase, dom::CartesianDomain2D, d::Real;
               fftshifted::Bool=false, coordmap=identity) -> ModalPhase

Convert a symbolic Zernike phase directly to a `ModalPhase` whose basis contains
exactly `ph`'s own modes (via [`zernike_basis`](@ref)) — no zero-padding, unlike
`ModalPhase(ph, basis::ZernikeBW)` which requires a pre-built basis containing
(at least) all of `ph`'s modes and pads the rest with zero coefficients.
"""
function ModalPhase(
    ph::SymbolicZernikePhase,
    dom::CartesianDomain2D,
    d::Real;
    fftshifted::Bool=false,
    coordmap=identity,
)
    basis = zernike_basis(ph, dom, d; fftshifted=fftshifted, coordmap=coordmap)
    return ModalPhase(copy(ph.coef), basis)
end
