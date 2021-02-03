
"""
    aperture(xrange, yrange, d, o=0)
    aperture(dom::CartesianDomain2D, d, o)

Create circular aperture in array `xrange × yrange` or given by `dom`ain with diameter `d` and center at `o`.

# Example
```jdoctests
julia> dom1= PhaseBases.CartesianDomain2D(-1:.2:1, -.9:.2:.9)
SampledDomains.CartesianDomain2D(-1.0:0.2:1.0, -0.9:0.2:0.9)

julia> ap, mask = PhaseBases.aperture(dom1, 1);

julia> ap
10×11 Array{Int64,2}:
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  1  0  0  0  0  0
 0  0  0  1  1  1  1  1  0  0  0
 0  0  0  1  1  1  1  1  0  0  0
 0  0  0  1  1  1  1  1  0  0  0
 0  0  0  1  1  1  1  1  0  0  0
 0  0  0  0  0  1  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0
```
"""
function aperture(xrange::AbstractRange, yrange::AbstractRange, d, o=(0,0))
ap = [ ((xc-o[1])^2 + (yc-o[2])^2) <= d^2 /4 ? 1 : 0 for yc ∈ yrange,  xc ∈ xrange]
# area = +(ap[:]...)
phmask = [ ((xc-o[1])^2 + (yc-o[2])^2) <= d^2 /4 ? 1 : NaN for yc ∈ yrange,  xc ∈ xrange]
return(ap, phmask)
end

aperture(dom::CartesianDomain2D, d, o=(0,0)) = aperture(dom.xrange, dom.yrange, d,o)