module PhaseBases
using  RecursiveArrayTools

export zernike, ZernikeBW

include("types.jl")


include("Zernike/Zernike.jl")

test(x) = x^3 



end

