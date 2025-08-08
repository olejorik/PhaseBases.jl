# it's supposed that basis is a vector of its elements
# Each element of the basis can be a multidimensional array ()
# Internally it can be implemented differently: as one multidimensional array or as vector of arrays.
# Externally it should make no difference via the following organization.
using RecursiveArrayTools
using SparseArrays

import Base: length, collect, promote_rule, convert, copy
import Base: -, +, *, show
import LinearAlgebra: norm, dot

export AbstractBasis, OrthogonalBasis, OrthonormalBasis, Phase, ModalPhase, ZonalPhase
export elements, aperture, aperturedelements, norms, dualelements, indexes, decompose, mask


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
function compose(b::AbstractBasis, ind::Vector, coef::Vector)
    return if length(ind) != length(coef)
        error("Coefficients and indexes are of different length")
    else
        comb(coef, elements(b, ind))
    end
end

function compose(b::AbstractBasis, coef::Vector)
    return if length(elements(b)) != length(coef)
        error("Coefficient vector does not match length of basis")
    else
        comb(coef, elements(b))
    end
end

@doc raw"""

    decompose(a, b::Basis)

Calculate coefficients of `a` in basis `b`.
"""
function decompose(a, b::AbstractBasis)
    if hasproperty(b, :dualelements)
        if hasproperty(b, :indexes)
            s = size(a)
            return [inner(a[CartesianIndex.(b.indexes)], f) for f in b.dualelements]
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

# Function below is for basis implemented as VectorOfArray
# Do we need the same for multidimensional array?
function comb(coef::Vector{T} where {T<:Number}, a::Union{VectorOfArray,AbstractVector})
    sum = similar(a[1])
    sum .= 0
    for i in 1:length(coef)
        sum += coef[i] * a[i]
    end
    return sum
end
comb(a::Union{VectorOfArray,AbstractVector}, coef::Vector{T} where {T<:Number}) =
    comb(coef, a)

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
inner(a, b, inds::Vector{T} where {T<:CartesianIndex}) = inner(a[inds], b[inds])


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
    coef::Vector{TC}
    basis::TB
    function ModalPhase{TC,TB}(
        coef::Vector{TC}, basis::TB
    ) where {TC<:Real,TB<:AbstractBasis}
        return if length(coef) != length(basis)
            error("Length of coefficent vector doesn't match the basis lenghth")
        else
            new(coef, basis)
        end
    end
end

function ModalPhase(coef::Vector{TC}, basis::TB) where {TC<:Real,TB<:AbstractBasis}
    return ModalPhase{TC,TB}(coef, basis)
end

ModalPhase(basis::AbstractBasis) = ModalPhase(zeros(Float64, length(basis)), basis)

function ModalPhase(
    ind::Vector{Int}, coef::Vector{TC}, basis::TB
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

norm(ph::ModalPhase) = sqrt(sum((ph.coef .* norms(ph.basis)) .^ 2))


copy(ph::ModalPhase) = ModalPhase(copy(ph.coef), ph.basis)

collect(ph::ModalPhase) = compose(ph.basis, ph.coef)
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
