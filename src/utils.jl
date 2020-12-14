
"""
aperture(xrange, yrange, d, o=0)
aperture(dom::CartesianDomain2D, d, o)

Create circular aperture in array `xrange × yrange` with diameter `d` and center at `o`.
"""
function aperture(xrange::AbstractRange, yrange::AbstractRange, d, o=(0,0))
ap = [ ((xc-o[1])^2 + (yc-o[2])^2) <= d^2 /4 ? 1 : 0 for yc ∈ yrange,  xc ∈ xrange]
# area = +(ap[:]...)
phmask = [ (xc^2 + yc^2) <= d^2 /4 ? 1 : NaN for yc ∈ yrange,  xc ∈ xrange]
return(ap, phmask)
end

aperture(dom::CartesianDomain2D, d, o=(0,0)) = aperture(dom.xrange, dom.yrange, d,o)