# it's supposed that basis is a vector of its elements
using RecursiveArrayTools

abstract type Basis end

abstract type OrthogonalBasis <: Basis end

abstract type OrthonormalBasis <: OrthogonalBasis end

"""

    elements(b) gives array of the basis functions.
"""
elements(b::Basis) = @view b.elements[:]

elements(b::Basis, ind::Vector) = @view b.elements[ind]

norms(b::Basis) = @view b.norms[:]

norms(b::Basis, ind::Vector) = @view b.norms[ind]



@doc raw"""
    compose(coef, ind, b) gives linear combination  math`\sum_i \lamda_i f_i`
"""
compose(coef::Vector, ind::Vector, b::Basis) = comb(coef, elements(b,ind))

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

#Function below is oversimplified, uses assumptions that b is 3d array
function comb(coef::Vector, a::Array)
    sum = similar(a[1]) .* 0
    for i = 1: length(coef)
        sum += coef[i] * a[i]
    end
    return sum
end


# a proper version of comb
"""
    inner(a, b)

Calculate inner product (tensor convolution) using the last index in `a` and first index in `b`
"""
function inner(a::Array, b::Array) end