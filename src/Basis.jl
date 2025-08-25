# General basis construction


struct Basis <: AbstractBasis
    elements::VectorOfArray
    dualelements::VectorOfArray
    ap::Array
    indexes::Array{Tuple}
    norms::Vector
    function Basis(elements, indexes; atol=0, rtol=0)
        ap = zeros(size(first(elements)))
        ap[CartesianIndex.(indexes)] .= 1
        elten = zeros(length(indexes), length(elements))
        for (i, e) in enumerate(elements)
            elten[:, i] = e[CartesianIndex.(indexes)]
        end
        # elten = reshape(Array(elements[CartesianIndex.(indexes)]), (:, length(elements)))
        if atol == 0 && rtol == 0
            invels = pinv(elten; atol=atol)
        else
            invels = pinv(elten; atol=atol, rtol=rtol)
        end
        dualelements = VectorOfArray(eachrow(invels))
        return new(
            elements, dualelements, ap, indexes, [sqrt.(inner(f, f)) for f in elements]
        )
    end
end

Basis(elements::Vector, indexes; kwargs...) =
    Basis(VectorOfArray(elements), indexes; kwargs...)

# We also introduce a basis with shifted origin
#= struct ShiftedBasis <: AbstractBasis
    elements::VectorOfArray
    origin::Array
    dualelements::VectorOfArray
    ap::Array
    indexes::Array{Tuple}
    norms::Vector
    function Basis(elements, origin, indexes; atol=0, rtol=0)
        ap = zeros(size(first(elements)))
        ap[indexes] .= 1
        elten = zeros(length(indexes), length(elements))
        for (i, e) in enumerate(elements)
            elten[:, i] = e[indexes]
        end
        # elten = reshape(Array(elements[CartesianIndex.(indexes)]), (:, length(elements)))
        if atol == 0 && rtol == 0
            invels = pinv(elten; atol=atol)
        else
            invels = pinv(elten; atol=atol, rtol=rtol)
        end
        dualelements = VectorOfArray(eachrow(invels))
        return new(
            elements,
            origin,
            dualelements,
            ap,
            indexes,
            [sqrt.(inner(f, f)) for f in elements],
        )
    end
end
 =#

struct ShiftedBasis <: AbstractBasis
    basis::Basis
    origin::Array
end

ShiftedBasis(elements, origin, indexes; kwargs...) =
    ShiftedBasis(Basis(elements, indexes; kwargs...), origin)

decompose(a, b::ShiftedBasis) = decompose(a .- b.origin, b.basis)
compose(b::ShiftedBasis, ind::Vector, coef::Vector) =
    compose(b.basis, ind, coef) .+ b.origin
compose(b::ShiftedBasis, coef::Vector) = compose(b.basis, coef) .+ b.origin

#TODO forward all other methods using https://github.com/curtd/ForwardMethods.jl
