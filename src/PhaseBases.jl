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
    compose!,
    decompose,
    decompose!,
    collect!,
    project,
    decompose_and_complement,
    ShiftedBasis
export zernike
export SymbolicZernikePhase,
    ZernikeOrdering, ZernikeNormalization,
    FringeOrdering, NollOrdering, OSAOrdering, MizerOrdering,
    BWNormalization, RMSNormalization,
    Fringe, Noll, OSA, Mizer, BornWolf, RMSNorm,
    reorder, to_nm, j_to_nm, nm_to_j, j_first

include("types.jl")
include("utils.jl")
include("Basis.jl")

include("Zernike/Zernike.jl")

end
