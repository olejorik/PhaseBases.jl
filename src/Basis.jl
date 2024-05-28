# General basis construction


struct Basis <: AbstractBasis
    elements::VectorOfArray
    dualelements::VectorOfArray
    ap::Array
    indexes::Array{Tuple}
    norms::Vector
    function Basis(elements, indexes)
        ap = ones(size(first(elements)))
        elten = reshape(Array(elements), (:, length(elements)))
        invels = pinv(elten)
        dualelements = VectorOfArray(eachrow(invels))
        return new(
            elements, dualelements, ap, indexes, [sqrt.(inner(f, f)) for f in elements]
        )
    end
end
