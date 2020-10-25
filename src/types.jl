# it's supposed that basis is a vector of its elements
# Each element of the basis can be a multidimensional array ()
# Internally it can be implemented differently: as one multidimensional array or as vector of arrays.
# Externally it should make no difference via the following organization.
using RecursiveArrayTools

abstract type Basis end

abstract type OrthogonalBasis <: Basis end

abstract type OrthonormalBasis <: OrthogonalBasis end

"""

    elements(b::Basis) gives vector of the basis functions.
"""
elements(b::Basis) =  b.elements[:]

elements(b::Basis, ind) =  b.elements[ind]

norms(b::Basis) =  b.norms[:]

norms(b::Basis, ind::Vector) =  b.norms[ind]



@doc raw"""
    compose(b, ind, coef) gives linear combination  math`\sum_i \lambda_i f_i`
"""
compose(b::Basis, ind::Vector, coef::Vector) = comb(coef, elements(b,ind))

@doc raw"""

    decompose(a::Array, b::Basis)

Calculate coefficients of `a` in basis `b`.
"""
function decompose(a::Array, b::Basis)
    error("There are no decomposition rules of basis $b") 
end

function decompose(a::Array, b::OrthonormalBasis)
    [ a[:] ⋅ f[:] for f in elements(b) ]
end

function decompose(a::Array, b::OrthogonalBasis)
    [ a[:] ⋅ f[:] / n for (f,n) in zip(elements(b),norms(b)) ]
end

# Function below is for basis implemented as VectorOfArray
# Do we need the same for multidimensional array?
function comb(coef::Vector, a::VectorOfArray)
    sum = similar(a[1]) .* 0
    for i = 1: length(coef)
        sum += coef[i] * a[i]
    end
    return sum
end


# a proper version of comb
# It can be realised as inner product in Mathematica
# Calculate inner product (tensor convolution) using the last index in `a` and first index in `b`
# Can be also related to scalar product in Banach space. Then inner product `(c, f)` of function `f` and
# constant `c` is equivalent to representing `c` as a constant function.
"""
    inner(a, b)

Calculate inner product of functions or array of functions.
"""
inner(a::Array{T,N}, b::Array{T,N}) where {T <: Number, N} = dot(a[:], b[:])

inner(a::Number, b::Number) = dot(a,b)
inner(c::Number, a::Array) = dot(c, sum(a))
inner(a::Array, c::Number) = dot(sum(a), c)

inner(a::Union{Vector, VectorOfArray}, b::Union{Vector, VectorOfArray}) = sum(inner(va,vb) for (va,vb) in zip(a, b))

innermatrix(b::Basis) = [inner(elements(b,i), b.ap .* elements(b,j)) for i in eachindex(elements(b)), j in eachindex(elements(b))]