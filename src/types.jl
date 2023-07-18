# it's supposed that basis is a vector of its elements
# Each element of the basis can be a multidimensional array ()
# Internally it can be implemented differently: as one multidimensional array or as vector of arrays.
# Externally it should make no difference via the following organization.
using RecursiveArrayTools

import Base: length, collect

export Basis, OrthogonalBasis, OrthonormalBasis, Phase, ModalPhase
export elements, aperture, aperturedelements, norms

"""
Abstract type representing any set of functions (`elements`) defined on some subset `ap` of a Cartesian domain.

    `elements(b::Basis)` gives vector of the basis functions.
    `aperture(b::Basis)` gives the  integer mask of non-zero elements of the array.
    `length(b::Basis)` gives the total number of elements in the basis.
"""
abstract type Basis end

abstract type OrthogonalBasis <: Basis end

abstract type OrthonormalBasis <: OrthogonalBasis end

elements(b::Basis) = b.elements
elements(b::Basis, ind) = b.elements[ind]

aperture(b::Basis) = b.ap

# Todo: correct below, returns array instead of VectorOfArray
aperturedelements(b::Basis) = (b.ap) .* b.elements
aperturedelements(b::Basis, ind) = (b.ap) .* b.elements[ind]

norms(b::Basis) = b.norms[:]

norms(b::Basis, ind::Vector) = b.norms[ind]

length(b::Basis) = length(elements(b))

@doc raw"""
    compose(b, ind, coef) gives linear combination  math`\sum_{i\in ind \lambda_i f_i`


    compose(b, coef) requires `coef` be the complete vector of the coefficients.
"""
function compose(b::Basis, ind::Vector, coef::Vector)
    return if length(ind) != length(coef)
        error("Coefficients and indexes are of different length")
    else
        comb(coef, elements(b, ind))
    end
end
function compose(b::Basis, coef::Vector)
    return if length(elements(b)) != length(coef)
        error("Coefficient vector does not match length of basis")
    else
        comb(coef, elements(b))
    end
end

@doc raw"""

    decompose(a::Array, b::Basis)

Calculate coefficients of `a` in basis `b`.
"""
function decompose(a::Array, b::Basis)
    return error("There are no decomposition rules of basis $typeof(b)")
end

function decompose(a::Array, b::OrthonormalBasis)
    return [inner(a, f, aperture(b)) for f in elements(b)]
end

function decompose(a::Array, b::OrthogonalBasis)
    return [inner(a, f, aperture(b)) / n^2 for (f, n) in zip(elements(b), norms(b))]
end

# Function below is for basis implemented as VectorOfArray
# Do we need the same for multidimensional array?
function comb(coef::Vector, a::VectorOfArray)
    sum = similar(a[1]) .* 0
    for i in 1:length(coef)
        sum += coef[i] * a[i]
    end
    return sum
end

# TODO rewrite so it works as tensor inner product if applied to two multidimensional arrays (or use Tensors.jl/ TensorOperations.jl?)
# a proper version of comb
# It can be realised as inner product in Mathematica
# Calculate inner product (tensor contraction) using the last index in `a` and first index in `b`
# Can be also related to scalar product in Banach space. Then inner product `(c, f)` of function `f` and
# constant `c` is defined as inner of `(f) ⋅ (c) = c f` .
"""
    inner(a, b)

Calculate inner product of functions or array of functions.
"""
inner(a::Array{T,N}, b::Array{T,N}) where {T<:Number,N} = integrate(a .* b)
inner(a::Array{T,N}, b::Array{T,N}, domain) where {T<:Number,N} = integrate(a .* b, domain)

inner(a::Number, b::Number) = dot(a, b)
# inner(c::Number, a::Array) = dot(c, sum(a))
# inner(a::Array, c::Number) = dot(sum(a), c)
inner(c::Number, a::Array) = conj(c) * a
inner(a::Array, c::Number) = conj(a) * c

function inner(a::Union{Vector,VectorOfArray}, b::Union{Vector,VectorOfArray})
    return sum(inner(va, vb) for (va, vb) in zip(a, b))
end

function innermatrix(b::Basis)
    return [
        inner(elements(b, i), b.ap .* elements(b, j)) for i in eachindex(elements(b)),
        j in eachindex(elements(b))
    ]
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
"""
abstract type Phase end
abstract type AbstractModalPhase <: Phase end
abstract type AbstractZonalPhase <: Phase end

"""
    phase = ModalPhase([coefficients,] basis)

Create phase represented by its `coefficients` in `basis`.
If no coefficents vector is given, use all zeros.

The representation is "lazy", so no data are copied. Use `collect` to obtain
an array representing the phase values.

Coefficients can be set through `phase.coef` field.

"""
struct ModalPhase{TC<:Real,TB<:Basis} <: AbstractModalPhase
    coef::Vector{TC}
    basis::TB
    function ModalPhase{TC,TB}(coef::Vector{TC}, basis::TB) where {TC<:Real,TB<:Basis}
        return if length(coef) != length(basis)
            error("Length of coefficent vector doesn't match the basis lenghth")
        else
            new(coef, basis)
        end
    end
end

function ModalPhase(coef::Vector{TC}, basis::TB) where {TC<:Real,TB<:Basis}
    return ModalPhase{TC,TB}(coef, basis)
end
ModalPhase(basis::Basis) = ModalPhase(zeros(Float64, length(basis)), basis)
function ModalPhase(
    ind::Vector{Int}, coef::Vector{TC}, basis::TB
) where {TC<:Real,TB<:Basis}
    c = zeros(TC, length(basis))
    c[ind] .= coef
    return ModalPhase{TC,TB}(c, basis)
end

collect(ph::ModalPhase) = compose(ph.basis, ph.coef)

import Base: -, +, *
-(ph::ModalPhase) = ModalPhase(-ph.coef, ph.basis)
*(c::Real, ph::ModalPhase) = ModalPhase(c * ph.coef, ph.basis)
*(ph::ModalPhase, c::Real) = *(c::Real, ph::ModalPhase)
function +(ph1::ModalPhase, ph2::ModalPhase)
    return if ph1.basis === ph2.basis
        ModalPhase(ph1.coef + ph2.coef, ph1.basis)
    else
        error("Cannot add modal phases in different bases")
    end
end

-(ph1::ModalPhase, ph2::ModalPhase) = ph1 + (-ph2)
