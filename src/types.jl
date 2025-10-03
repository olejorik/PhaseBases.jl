# it's supposed that basis is a vector of its elements
# Each element of the basis can be a multidimensional array ()
# Internally it can be implemented differently: as one multidimensional array or as vector of arrays.
# Externally it should make no difference via the following organization.
using RecursiveArrayTools
using SparseArrays

import Base: length, collect, promote_rule, convert, copy
import Base: -, +, *, show
import LinearAlgebra: norm, dot

export AbstractBasis, OrthogonalBasis, OrthonormalBasis, Phase, ModalPhase, ZonalPhase, PixelBasis
export elements, aperture, aperturedelements, norms, dualelements, indexes, decompose, mask, coefficients, coefficients!


"""
Abstract type representing any set of functions (`elements`) defined on some subset `ap` of a Cartesian domain.

    `elements(b::Basis)` gives vector of the basis functions.
    `norms(b::Basis)` gives vector of the norms of the basis functions.
    `aperture(b::Basis)` gives the  integer mask of non-zero elements of the array.
    `length(b::Basis)` gives the total number of elements in the basis.
"""
abstract type AbstractBasis end

abstract type OrthogonalBasis <: AbstractBasis end

abstract type OrthonormalBasis <: OrthogonalBasis end

"""
    `elements(b::Basis [, ind])`

Gives vector of the basis functions (all of them or specified by the vector of integer indexes `ind`).
"""
elements(b::AbstractBasis) = b.elements
elements(b::AbstractBasis, ind) = b.elements[ind]

aperture(b::AbstractBasis) = b.ap

## aperture mask interface
"""
    mask(b::AbstractBasis)

Return the plotting mask array for basis `b`.
"""
function mask(b::AbstractBasis)
    if hasproperty(b, :mask)
        return getfield(b, :mask)
    else
        error("Basis of type $(typeof(b)) does not have a `mask` field")
    end
end

# aperturedelements: mask each basis function by aperture or mask
"""
    aperturedelements(b::AbstractBasis; use_mask=false)
    aperturedelements(b::AbstractBasis, ind::Integer; use_mask=false)
    aperturedelements(b::AbstractBasis, inds::AbstractVector{Int}; use_mask=false)

Return the basis elements masked by `aperture(b)` (default) or by `mask(b)` if `use_mask=true`.
With no index returns a vector; with integer returns a single array; with vector returns a vector of arrays.
"""
function aperturedelements(b::AbstractBasis; use_mask=false)
    mat = use_mask ? mask(b) : aperture(b)
    return [mat .* f for f in elements(b)]
end
function aperturedelements(b::AbstractBasis, ind::Integer; use_mask=false)
    mat = use_mask ? mask(b) : aperture(b)
    return mat .* elements(b, ind)
end
function aperturedelements(b::AbstractBasis, inds::AbstractVector{Int}; use_mask=false)
    mat = use_mask ? mask(b) : aperture(b)
    return [mat .* f for f in elements(b, inds)]
end

norms(b::AbstractBasis) = b.norms[:]

norms(b::AbstractBasis, ind::Vector) = b.norms[ind]

length(b::AbstractBasis) = length(elements(b))

# indexing shortcut: allow b[i] or b[inds]
import Base: getindex
getindex(b::AbstractBasis, i::Integer) = elements(b, i)
getindex(b::AbstractBasis, inds::AbstractVector{Int}) = elements(b, inds)

"""
    dualelements(b::AbstractBasis)

Return the dual (pseudo-inverse) elements of basis `b`.
"""
function dualelements(b::AbstractBasis)
    if hasproperty(b, :dualelements)
        return getfield(b, :dualelements)
    else
        error("Basis of type $(typeof(b)) does not support dualelements")
    end
end

"""
    indexes(b::AbstractBasis)

Return the pixel indices where basis `b` is defined (only for sparse or indexed bases).
"""
function indexes(b::AbstractBasis)
    if hasproperty(b, :indexes)
        return getfield(b, :indexes)
    else
        error("Basis of type $(typeof(b)) does not have an `indexes` field")
    end
end

@doc raw"""
    compose(b, ind, coef) gives linear combination  math`\sum_{i\in ind \lambda_i f_i`


    compose(b, coef) requires `coef` be the complete vector of the coefficients.
"""
function compose(b::AbstractBasis, ind::AbstractVector, coef::AbstractVector)
    return if length(ind) != length(coef)
        error("Coefficients and indexes are of different length")
    else
        comb(coef, elements(b, ind))
    end
end

function compose(b::AbstractBasis, coef::AbstractVector)
    return if length(elements(b)) != length(coef)
        error("Coefficient vector does not match length of basis")
    else
        comb(coef, elements(b))
    end
end

compose!(target, b::AbstractBasis, coef::AbstractVector) =
    (size(target) == size(aperture(b)) ||
        throw(ArgumentError("Target array size does not match basis aperture size"));
    comb!(target, coef, elements(b)))

compose!(target, b::AbstractBasis, ind::AbstractVector, coef::AbstractVector) =
    (size(target) == size(aperture(b)) ||
        throw(ArgumentError("Target array size does not match basis aperture size"));
    comb!(target, coef, elements(b, ind)))



"""
    decompose(a, b::Basis)

Calculate coefficients of `a` in basis `b`.
"""
function decompose(a, b::AbstractBasis)
    if hasproperty(b, :dualelements)
        if hasproperty(b, :indexes)
            return [inner(a, f, b.indexes) for f in b.dualelements]
        else #assume the indexes is the full array
            return [inner(a, f, aperture(b)) for f in b.dualelements]
        end
    else
        return error("There are no decomposition rules of basis $typeof(b)")
    end
end

function decompose(a, b::OrthonormalBasis)
    return [inner(a, f, aperture(b)) for f in elements(b)]
end

function decompose(a, b::OrthogonalBasis)
    return [inner(a, f, aperture(b)) / n^2 for (f, n) in zip(elements(b), norms(b))]
end

"""
    decompose!(coeffs, a, b::AbstractBasis)

Non-allocating version of `decompose`. Write coefficients of `a` in basis `b` into pre-allocated `coeffs` vector.
"""
function decompose!(coeffs::AbstractVector, a, b::AbstractBasis)
    length(coeffs) == length(b) ||
        throw(ArgumentError("Coefficient vector length $(length(coeffs)) doesn't match basis length $(length(b))"))

    if hasproperty(b, :dualelements)
        if hasproperty(b, :indexes)
            @inbounds for (i, f) in enumerate(b.dualelements)
                coeffs[i] = inner_indexed(a, f, b.indexes)
            end
        else #assume the indexes is the full array
            @inbounds for (i, f) in enumerate(b.dualelements)
                coeffs[i] = inner_masked(a, f, aperture(b))
            end
        end
    else
        error("There are no decomposition rules of basis $(typeof(b))")
    end
    return coeffs
end

function decompose!(coeffs::AbstractVector, a, b::OrthonormalBasis)
    length(coeffs) == length(b) ||
        throw(ArgumentError("Coefficient vector length $(length(coeffs)) doesn't match basis length $(length(b))"))

    @inbounds for (i, f) in enumerate(elements(b))
        coeffs[i] = inner_masked(a, f, aperture(b))
    end
    return coeffs
end

function decompose!(coeffs::AbstractVector, a, b::OrthogonalBasis)
    length(coeffs) == length(b) ||
        throw(ArgumentError("Coefficient vector length $(length(coeffs)) doesn't match basis length $(length(b))"))

    basis_norms = norms(b)
    @inbounds for (i, f) in enumerate(elements(b))
        coeffs[i] = inner_masked(a, f, aperture(b)) / basis_norms[i]^2
    end
    return coeffs
end

project(a, b::AbstractBasis) = compose(b, decompose(a, b))

decompose_and_complement(a, b::AbstractBasis) = (decompose(a, b), a .- project(a, b))

# These functions are needed for calculation of the derivatives by direction

"""
    allinners!(coeffs, a, b::AbstractBasis)

Write results of inner product of `a` with all elements of basis `b` into pre-allocated `coeffs` vector.
"""
function allinners!(coeffs::AbstractVector, a, b::AbstractBasis)
    length(coeffs) == length(b) ||
        throw(ArgumentError("Coefficient vector length $(length(coeffs)) doesn't match basis length $(length(b))"))

    if hasproperty(b, :elements)
        if hasproperty(b, :indexes)
            @inbounds for (i, f) in enumerate(b.elements)
                coeffs[i] = inner_indexed(a, f, b.indexes)
            end
        else #assume the indexes is the full array
            @inbounds for (i, f) in enumerate(b.elements)
                coeffs[i] = inner_masked(a, f, aperture(b))
            end
        end
    else
        error("There are no decomposition rules of basis $(typeof(b))")
    end
    return coeffs
end

function allinners!(coeffs::AbstractVector, a, b::OrthonormalBasis)
    length(coeffs) == length(b) ||
        throw(ArgumentError("Coefficient vector length $(length(coeffs)) doesn't match basis length $(length(b))"))

    @inbounds for (i, f) in enumerate(elements(b))
        coeffs[i] = inner_masked(a, f, aperture(b))
    end
    return coeffs
end

function allinners!(coeffs::AbstractVector, a, b::OrthogonalBasis)
    length(coeffs) == length(b) ||
        throw(ArgumentError("Coefficient vector length $(length(coeffs)) doesn't match basis length $(length(b))"))

    @inbounds for (i, f) in enumerate(elements(b))
        coeffs[i] = inner_masked(a, f, aperture(b)) 
    end
    return coeffs
end


# Function below is for basis implemented as VectorOfArray
# Do we need the same for multidimensional array?
function comb!(target, coef::AbstractVector{T} where {T<:Number}, a::Union{VectorOfArray,AbstractVector})
    size(target) == size(a[1]) ||
        throw(ArgumentError("Target array size does not match basis element size"))
    target .= 0
    for i in 1:length(coef)
        target .+= coef[i] * a[i]
    end
    return target
end

function comb(coef::AbstractVector{T} where {T<:Number}, a::Union{VectorOfArray,AbstractVector})
    sum = similar(a[1])
    # sum .= 0
    # for i in 1:length(coef)
    #     sum += coef[i] * a[i]
    # end
    # return sum
    comb!(sum, coef, a)
    return sum
end
comb(a::Union{VectorOfArray,AbstractVector}, coef::AbstractVector{T} where {T<:Number}) =
    comb(coef, a)

comb!(target, a::Union{VectorOfArray,AbstractVector}, coef::AbstractVector{T} where {T<:Number}) =
    comb!(target, coef, a)

# TODO rewrite so it works as tensor inner product if applied to two multidimensional arrays (or use Tensors.jl/ TensorOperations.jl?)
# a proper version of comb
# It can be realised as inner product in Mathematica
# Calculate inner product (tensor contraction) using the last index in `a` and first index in `b`
# Can be also related to scalar product in Banach space. Then inner product `(c, f)` of function `f` and
# constant `c` is defined as inner of `(f) ⋅ (c) = c f` .
"""
    inner(a, b)

Calculate inner product of sampled functions or array of functions.
"""
inner(a, b) = _inner(a, b)
_inner(a, b) = dot(a, b)

# The def above is the same what is below but simpler and is general
# Expressed through the unsafe `_inner` function
# inner(a::Array{T,N}, b::Array{T,N}) where {T<:Number,N} = integrate(a .* b)
# inner(a::Number, b::Number) = dot(a, b)
# inner(c::Number, a::Array) = conj(c) * a
# inner(a::Array, c::Number) = conj(a) * c
#

inner(a, b, weights) = inner(a .* weights, b)
## Safely gather values at Vector{CartesianIndex} without relying on A[inds]
inner(a, b, inds::Vector{T}) where {T<:CartesianIndex} = inner(a[inds], b[inds])
## If `a` is already a vector over the aperture points, only index `b`
inner(a::AbstractVector, b, inds::Vector{T}) where {T<:CartesianIndex} = inner(a, b[inds])
## If `b` is already a vector over the aperture points, only index `a`
inner(a::AbstractArray, b::AbstractVector, inds::Vector{T}) where {T<:CartesianIndex} =
    dot(a[inds], b)
## Both arrays with Cartesian indices -> gather both and dot
inner(a::AbstractArray, b::AbstractArray, inds::Vector{T}) where {T<:CartesianIndex} =
    dot(a[inds], b[inds])

# ===============================================================================
# Allocation-free inner product functions
# ===============================================================================

"""
    inner_masked(a, b, mask)

Allocation-free inner product of `a` and `b` with boolean `mask`.
Only computes dot product at locations where mask is true.
"""
function inner_masked(a::AbstractArray, b::AbstractArray, mask::AbstractArray{Bool})
    size(a) == size(b) == size(mask) ||
        throw(ArgumentError("Arrays must have same size"))

    result = zero(promote_type(eltype(a), eltype(b)))
    @inbounds for i in eachindex(a, b, mask)
        if mask[i]
            result += conj(a[i]) * b[i]
        end
    end
    return result
end

"""
    inner_weighted(a, b, weights)

Allocation-free weighted inner product. Equivalent to `inner(a .* weights, b)` 
but without temporary array allocation.
"""
function inner_weighted(a::AbstractArray, b::AbstractArray, weights::AbstractArray)
    size(a) == size(b) == size(weights) ||
        throw(ArgumentError("Arrays must have same size"))

    result = zero(promote_type(eltype(a), eltype(b), eltype(weights)))
    @inbounds for i in eachindex(a, b, weights)
        result += conj(a[i]) * weights[i] * b[i]
    end
    return result
end

"""
    inner_indexed(a, b, inds::Vector{CartesianIndex})

Allocation-free inner product at specified indices. Equivalent to `dot(a[inds], b[inds])`
but without creating temporary arrays.
"""
function inner_indexed(a::AbstractArray, b::AbstractArray, inds::Vector{CartesianIndex{N}}) where {N}
    result = zero(promote_type(eltype(a), eltype(b)))
    @inbounds for idx in inds
        result += conj(a[idx]) * b[idx]
    end
    return result
end

"""
    inner_indexed(a::AbstractVector, b, inds::Vector{CartesianIndex})

Allocation-free inner product when `a` is already a vector over aperture points.
"""
function inner_indexed(a::AbstractVector, b::AbstractArray, inds::Vector{CartesianIndex{N}}) where {N}
    length(a) == length(inds) ||
        throw(ArgumentError("Vector length $(length(a)) doesn't match index length $(length(inds))"))

    result = zero(promote_type(eltype(a), eltype(b)))
    @inbounds for (i, idx) in enumerate(inds)
        result += conj(a[i]) * b[idx]
    end
    return result
end

"""
    inner_indexed(a::AbstractArray, b::AbstractVector, inds::Vector{CartesianIndex})

Allocation-free inner product when `b` is already a vector over aperture points.
"""
function inner_indexed(a::AbstractArray, b::AbstractVector, inds::Vector{CartesianIndex{N}}) where {N}
    length(b) == length(inds) ||
        throw(ArgumentError("Vector length $(length(b)) doesn't match index length $(length(inds))"))

    result = zero(promote_type(eltype(a), eltype(b)))
    @inbounds for (i, idx) in enumerate(inds)
        result += conj(a[idx]) * b[i]
    end
    return result
end

# ===============================================================================
# Enhanced decompose! functions using allocation-free inner products
# ===============================================================================

"""
    decompose_masked!(coeffs, a, b::OrthonormalBasis, mask::AbstractArray{Bool})

Non-allocating decompose using boolean mask instead of aperture multiplication.
"""
function decompose_masked!(coeffs::AbstractVector, a, b::OrthonormalBasis, mask::AbstractArray{Bool})
    length(coeffs) == length(b) ||
        throw(ArgumentError("Coefficient vector length $(length(coeffs)) doesn't match basis length $(length(b))"))

    @inbounds for (i, f) in enumerate(elements(b))
        coeffs[i] = inner_masked(a, f, mask)
    end
    return coeffs
end


function inner(a::AbstractSparseArray, b::AbstractSparseArray)
    (a.colptr == b.colptr && a.rowval == b.rowval) || return dot(a, b)
    return _inner(nonzeros(a), nonzeros(b))
end

"""
    _inner(a::AbstractSparseArray, b::AbstractSparseArray)

Calculate inner product without checking the position of nonzero elements.
"""
function _inner(a::AbstractSparseArray, b::AbstractSparseArray)
    return dot(nonzeros(a), nonzeros(b))
end

function inner(a::Union{Vector,VectorOfArray}, b::Union{Vector,VectorOfArray})
    # TODO rework with preallocated
    return sum(a[i] * b[i] for i in eachindex(a))
end

# function innermatrix(b::Basis)
#     return [
#         _inner(elements(b, i), b.ap .* elements(b, j)) for i in eachindex(elements(b)),
#         j in eachindex(elements(b))
#     ]
# end

"""
    innermatrix(
    a::Union{Vector,VectorOfArray}, b::Union{Vector,VectorOfArray}, weight=1
)

@doctest
```julia
julia> PhaseBases.innermatrix([10, 100], [2,3,4])
3×2 Matrix{Int64}:
 20  200
 30  300
 40  400
```
"""
function innermatrix(
    a::Union{Vector,VectorOfArray}, b::Union{Vector,VectorOfArray}, weight=1
)
    return [_inner(b[i], weight .* a[j]) for i in eachindex(b), j in eachindex(a)]
end

innermatrix(b::AbstractBasis) = innermatrix(elements(b), elements(b), aperture(b))

function dualbasis(b::AbstractBasis)
    ata = innermatrix(b)
    return dualel = [inner(elements(b), (inv(ata))[:, i]) for i in 1:length(elements(b))]
end

"""
    integrate(a)

calculate numerical integral of a
"""
# rectangle rule
integrate(a::Array{T,N}) where {T<:Number,N} = sum(a[:])
integrate(a::Array{T,N}, domain::Array{T,N}) where {T<:Number,N} = sum(a[:] .* domain[:])

# TODO simpson's rule on unit disk

"""
`Phase` is an abstract type for various phase descriptions and models.

The actual phase values are kept in `coef` field.
"""
abstract type Phase end
abstract type AbstractModalPhase <: Phase end
abstract type AbstractZonalPhase <: Phase end

# Interfaces

coefficients(ph::Phase) = error("No `coefficients` method is defined for $(typeof(ph))")
coefficients!(ph::Phase, coef) =
    error("No `coefficients!` method is defined for $(typeof(ph))")
# Add getting and setting index of the coefficient

# In place scaling coefficients through indexing
# scalecoefficients!(ph::Phase, c::Real) =

-(ph::Phase) = (ph1 = copy(ph); (ph1.coef .*= -1); ph1)
*(c::Real, ph::Phase) = (ph1 = copy(ph); (ph1.coef .*= c); ph1)
*(ph::Phase, c::Real) = *(c::Real, ph::Phase)
+(x::Phase, y::Phase) = +(promote(x, y)...)
-(x::Phase, y::Phase) = -(promote(x, y)...)
struct ZonalPhase <: AbstractZonalPhase
    coef::Array{Float64,2}
end

coefficients(ph::ZonalPhase) = ph.coef
function coefficients!(ph::ZonalPhase, coef)
    ph.coef .= coef
    return ph
end

getindex(ph::ZonalPhase, ind...) = getindex(coefficients(ph), ind...)

copy(ph::ZonalPhase) = ZonalPhase(copy(ph.coef))
collect(ph::ZonalPhase) = ph.coef

+(x::ZonalPhase, y::ZonalPhase) = ZonalPhase(x.coef + y.coef)
-(x::ZonalPhase, y::ZonalPhase) = ZonalPhase(x.coef - y.coef)

"""
    phase = ModalPhase([coefficients,] basis)

Create phase represented by its `coefficients` in `basis`.
If no coefficents vector is given, use all zeros.

The representation is "lazy", so no data are copied. Use `collect` to obtain
an array representing the phase values.

Coefficients can be set through `phase.coef` field.

"""
struct ModalPhase{TC<:Real,TB<:AbstractBasis} <: AbstractModalPhase
    coef::AbstractVector{TC}
    basis::TB
    function ModalPhase{TC,TB}(
        coef::AbstractVector{TC}, basis::TB
    ) where {TC<:Real,TB<:AbstractBasis}
        return if length(coef) != length(basis)
            error("Length of coefficent vector doesn't match the basis lenghth")
        else
            new(coef, basis)
        end
    end
end
function ModalPhase(coef::AbstractVector{TC}, basis::TB) where {TC<:Real,TB<:AbstractBasis}
    return ModalPhase{TC,TB}(coef, basis)
end

ModalPhase(basis::AbstractBasis) = ModalPhase(zeros(Float64, length(basis)), basis)

function ModalPhase(
    ind::AbstractVector{Int}, coef::AbstractVector{TC}, basis::TB
) where {TC<:Real,TB<:AbstractBasis}
    c = zeros(TC, length(basis))
    c[ind] .= coef
    return ModalPhase{TC,TB}(c, basis)
end

coefficients(ph::ModalPhase) = ph.coef
function coefficients!(ph::ModalPhase, coef)
    ph.coef .= coef
    return ph
end

getindex(ph::ModalPhase, ind...) = sum(
    coefficients(ph)[i] * getindex(elements(ph.basis)[i], ind...) for
    i in eachindex(elements(ph.basis))
)



norm(ph::ModalPhase) = sqrt(sum((ph.coef .* norms(ph.basis)) .^ 2)) # TODO change, this is valid only for orthogonal basis


copy(ph::ModalPhase) = ModalPhase(copy(ph.coef), ph.basis)

collect(ph::ModalPhase) = compose(ph.basis, ph.coef)

collect!(target, ph::ModalPhase) = compose!(target, ph.basis, ph.coef)



convert(::Type{ZonalPhase}, ph::ModalPhase) = ZonalPhase(collect(ph))
ZonalPhase(ph::ModalPhase) = ZonalPhase(collect(ph))
promote_rule(::Type{ModalPhase{T,B}}, ::Type{ZonalPhase}) where {T,B} = ZonalPhase

function +(ph1::ModalPhase, ph2::ModalPhase)
    return if ph1.basis === ph2.basis
        ModalPhase(ph1.coef + ph2.coef, ph1.basis)
    else
        error("Cannot add modal phases in different bases")
    end
end

# -(ph1::Phase, ph2::Phase) = ph1 + (-ph2)

# pretty printing
function show(io::IO, x::ModalPhase{TC,TB}) where {TC,TB}
    return println(io, "Modal phase in $TB basis with coefficients $(x.coef)")
end

# ===============================================================================
# PixelBasis - Computational pixel-based basis functions
# ===============================================================================

"""
    PixelBasis{T,N} <: OrthonormalBasis

A computational pixel-based basis where each basis function is a unit function
at a specific pixel location within the aperture. This basis is orthonormal and
highly efficient for pixel-based phase representations.

Each basis function is conceptually an array of zeros with a single unit element
at one of the aperture pixel locations, but these are computed on-demand rather
than stored to save memory.

# Fields
- `ap::Array{T,N}`: The aperture mask defining valid pixel locations
- `indexes::Vector{CartesianIndex{N}}`: Valid pixel locations within the aperture
- `size::NTuple{N,Int}`: Dimensions of the aperture array

# Example
```julia
# Create circular aperture
aperture = create_circular_aperture(64, 64)
basis = PixelBasis(aperture)

# Use with ModalPhase
coeffs = randn(length(basis))
phase = ModalPhase(coeffs, basis)
phase_array = collect(phase)  # Very efficient conversion
```
"""
struct PixelBasis{T,N} <: OrthonormalBasis
    ap::Array{T,N}
    indexes::Vector{CartesianIndex{N}}
    size::NTuple{N,Int}
end

# ===============================================================================
# Phase 1: Constructors
# ===============================================================================

"""
    PixelBasis(aperture::Array{T,N}) where {T,N}

Create a PixelBasis from an aperture mask. Non-zero elements of the aperture
define the valid pixel locations for basis functions.
"""
function PixelBasis(aperture::Array{T,N}) where {T,N}
    indexes = findall(!iszero, aperture)
    return PixelBasis{T,N}(aperture, indexes, size(aperture))
end

"""
    PixelBasis(indexes::Vector{CartesianIndex{N}}, size::NTuple{N,Int}) where {N}

Create a PixelBasis from explicit pixel indexes and array size.
"""
function PixelBasis(indexes::Vector{CartesianIndex{N}}, size::NTuple{N,Int}) where {N}
    ap = zeros(Bool, size)
    ap[indexes] .= true
    return PixelBasis{Bool,N}(ap, indexes, size)
end

"""
    PixelBasis(mask::BitArray{N}) where {N}

Create a PixelBasis from a boolean mask.
"""
function PixelBasis(mask::BitArray{N}) where {N}
    indexes = findall(mask)
    return PixelBasis{Bool,N}(mask, indexes, size(mask))
end

# ===============================================================================
# Phase 3: Core Interface Methods (Optimized)
# ===============================================================================

"""
    compose(b::PixelBasis{T,N}, coef::Vector) where {T,N}

Efficiently compose pixel coefficients into a full array. This is the main
advantage of PixelBasis - O(N_pixels) composition via direct indexing.
"""
function compose(b::PixelBasis{T,N}, coef::AbstractVector) where {T,N}
    length(coef) == length(b.indexes) ||
        throw(ArgumentError("Coefficient length $(length(coef)) doesn't match basis length $(length(b.indexes))"))

    result = zeros(eltype(coef), b.size)
    result[b.indexes] .= coef
    return result
end

function compose!(target, b::PixelBasis{T,N}, coef::AbstractVector) where {T,N}
    size(target) == b.size ||
        throw(ArgumentError("Target array size $(size(target)) doesn't match basis size $(b.size)"))
    length(coef) == length(b.indexes) ||
        throw(ArgumentError("Coefficient length $(length(coef)) doesn't match basis length $(length(b.indexes))"))

    target .= 0
    target[b.indexes] .= coef
    return target
end

"""
    decompose(array::Array{T,N}, b::PixelBasis{T,N}) where {T,N}

Efficiently extract pixel coefficients from a full array. O(N_pixels) operation
via direct indexing.
"""
function decompose(array::Array{S,N}, b::PixelBasis{T,N}) where {S,T,N}
    size(array) == b.size ||
        throw(ArgumentError("Array size $(size(array)) doesn't match basis size $(b.size)"))

    return array[b.indexes]
end

"""
    decompose!(coeffs, array::Array{T,N}, b::PixelBasis{T,N}) where {T,N}

Non-allocating version of PixelBasis decompose. Write pixel coefficients from `array` 
into pre-allocated `coeffs` vector. O(N_pixels) operation via direct indexing.
"""
function decompose!(coeffs::AbstractVector, array::Array{S,N}, b::PixelBasis{T,N}) where {S,T,N}
    size(array) == b.size ||
        throw(ArgumentError("Array size $(size(array)) doesn't match basis size $(b.size)"))
    length(coeffs) == length(b.indexes) ||
        throw(ArgumentError("Coefficient vector length $(length(coeffs)) doesn't match basis length $(length(b.indexes))"))

    @inbounds for (i, idx) in enumerate(b.indexes)
        coeffs[i] = array[idx]
    end
    return coeffs
end

allinners!(coeffs::AbstractVector, array::Array{S,N}, b::PixelBasis{T,N}) where {S,T,N} = decompose!(coeffs, array, b)

# Required AbstractBasis interface methods
aperture(b::PixelBasis) = b.ap
Base.length(b::PixelBasis) = length(b.indexes)
norms(b::PixelBasis{T}) where {T} = ones(real(T), length(b.indexes))  # Orthonormal
indexes(b::PixelBasis) = b.indexes

# ===============================================================================
# Elements interface - computational implementation
# ===============================================================================

"""
    elements(b::PixelBasis)

Return a lazy generator that produces all basis functions on-demand. This is memory
efficient as it doesn't store all basis functions simultaneously.
"""
elements(b::PixelBasis) = (
    let
        result = zeros(eltype(b.ap), b.size)
        result[b.indexes[i]] = one(eltype(b.ap))
        result
    end for i in 1:length(b)
)

"""
    elements(b::PixelBasis, i::Int)

Generate the i-th basis function as an array. The basis function is a unit function
(array of zeros with a single 1) at the i-th pixel location within the aperture.
"""
function elements(b::PixelBasis{T,N}, i::Int) where {T,N}
    1 ≤ i ≤ length(b.indexes) || throw(BoundsError(b.indexes, i))

    result = zeros(T, b.size)
    result[b.indexes[i]] = one(T)
    return result
end

"""
    elements(b::PixelBasis, inds)

Generate multiple basis functions specified by indices `inds`.
"""
elements(b::PixelBasis, inds) = [elements(b, i) for i in inds]

# Pretty printing
function show(io::IO, x::PixelBasis{T,N}) where {T,N}
    return println(io, "PixelBasis{$T,$N} with $(length(x)) pixels in $(x.size) array")
end
