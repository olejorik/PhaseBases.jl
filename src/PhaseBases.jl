module PhaseBases
using RecursiveArrayTools
using LinearAlgebra
using SampledDomains: CartesianDomain2D
using FFTW

export Basis,
    ZernikeBW,
    ZernikeBWSparse,
    PixelBasis,
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
