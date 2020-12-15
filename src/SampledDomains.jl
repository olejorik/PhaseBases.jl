module SampledDomains

export CartesianDomain2D, make_centered_domain2D

#TODO rewrite as parametric type
struct CartesianDomain2D
    xrange::AbstractRange
    yrange::AbstractRange
end

function make_centered_domain2D(xlength, ylength, pixelsize)
    xrange = ((1:xlength) .- (1+xlength)/2) .* pixelsize
    yrange = ((1:ylength) .- (1+ylength)/2) .* pixelsize
    CartesianDomain2D(xrange, yrange)
end


end
