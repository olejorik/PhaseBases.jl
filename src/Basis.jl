# General basis construction


struct Basis <: AbstractBasis
    elements::VectorOfArray
    dualelements::VectorOfArray
    ap::Array
    indexes::Array{Tuple}
    norms::Vector
    function Basis(elements, indexes; atol=0, rtol=0)
        ap = ones(size(first(elements)))
        # ap[indexes] .= 1
        elten = reshape(Array(elements), (:, length(elements)))
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
