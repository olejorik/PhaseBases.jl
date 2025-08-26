module PhaseBases
using RecursiveArrayTools
using LinearAlgebra
using SampledDomains: CartesianDomain2D

export Basis,
    ZernikeBW,
    ZernikeBWSparse,
    compose,
    decompose,
    project,
    decompose_and_complement,
    ShiftedBasis
export zernike

include("types.jl")
include("utils.jl")
include("Basis.jl")

include("Zernike/Zernike.jl")

end
