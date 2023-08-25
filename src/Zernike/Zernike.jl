@doc raw"""
    zernike(x, y, maxord)

Calculate  unit-normalized Zernike polynomials and their x,y derivatives

Translation of the code from
[1] T. B. Andersen, “Efficient and robust recurrence relations for the Zernike circle polynomials
and their derivatives in Cartesian coordinates,” _Opt. Express_, vol. 26, no. 15, p. 18878, Jul. 2018.

     Numbering scheme:
     Within a radial order, sine terms come first
             ...
           sin((n-2m)*theta)   for m = 0,..., [(n+1)/2]-1
             ...
              1                for n even, m = n/2
             ...
           cos((n-2m)*theta)   for m = [n/2]+1,...,n
             ...

     INPUT:
     x, y normalized (x,y) coordinates in unit circle
     MaxOrd: Maximum Zernike radial order

     OUTPUT:
     dUdx[...]   array to receive each derivative dU/dx at (x,y)
     dUdy[...]   array to receive each derivative dU/dy at (x,y)

"""
function zernike(x, y, maxord::Int)
    ztotnumber = Int((maxord + 2) * (maxord + 1) / 2)
    # print(ztotnumber)
    zern = zeros(ztotnumber)
    dudx = zeros(ztotnumber)
    dudy = zeros(ztotnumber)

    zern[1] = 1
    dudx[1] = 0
    dudy[1] = 0

    zern[2] = y
    zern[3] = x
    dudx[2] = 0
    dudx[3] = 1
    dudy[2] = 1
    dudy[3] = 0

    kndx = 1                # index for term from 2 orders down
    jbeg = 2                # start index for current radial order
    jend = 3                # end index for current radial order
    jndx = 3                # running index for current Zern
    even = -1

    #  Outer loop in radial order index
    for nn in 2:maxord
        even = -even          # parity of radial index
        jndx1 = jbeg           # index for 1st ascending series in x
        jndx2 = jend           # index for 1st descending series in y
        jndx11 = jndx1 - 1      # index for 2nd ascending series in x
        jndx21 = jndx2 + 1      # index for 2nd descending series in y
        jbeg = jend + 1       # end of previous radial order +1
        nn2 = nn / 2
        nn1 = (nn - 1) / 2

        # Inner loop in azimuthal index
        for mm in 0:nn
            jndx += 1                  # increment running index for current Zern

            if mm == 0
                zern[jndx] = x * zern[jndx1] + y * zern[jndx2]
                dudx[jndx] = zern[jndx1] * nn
                dudy[jndx] = zern[jndx2] * nn

            elseif mm == nn
                zern[jndx] = x * zern[jndx11] - y * zern[jndx21]
                dudx[jndx] = zern[jndx11] * nn
                dudy[jndx] = -zern[jndx21] * nn

            elseif even > 0 && mm == nn2
                zern[jndx] = 2 * (x * zern[jndx1] + y * zern[jndx2]) - zern[kndx]
                dudx[jndx] = 2 * nn * zern[jndx1] + dudx[kndx]
                dudy[jndx] = 2 * nn * zern[jndx2] + dudy[kndx]
                kndx += 1

            elseif even < 0 && mm == nn1
                qval = zern[jndx2] - zern[jndx21]
                zern[jndx] = x * zern[jndx11] + y * qval - zern[kndx]
                dudx[jndx] = zern[jndx11] * nn + dudx[kndx]
                dudy[jndx] = qval * nn + dudy[kndx]
                kndx += 1

            elseif even < 0 && mm == nn1 + 1
                pval = zern[jndx1] + zern[jndx11]
                zern[jndx] = x * pval + y * zern[jndx2] - zern[kndx]
                dudx[jndx] = pval * nn + dudx[kndx]
                dudy[jndx] = zern[jndx2] * nn + dudy[kndx]
                kndx += 1

            else
                pval = zern[jndx1] + zern[jndx11]
                qval = zern[jndx2] - zern[jndx21]
                zern[jndx] = x * pval + y * qval - zern[kndx]
                dudx[jndx] = pval * nn + dudx[kndx]
                dudy[jndx] = qval * nn + dudy[kndx]
                kndx += 1
            end

            jndx11 = jndx1                   # update indices
            jndx1 += 1
            jndx21 = jndx2
            jndx2 += -1
        end

        jend = jndx
    end

    return (z=zern, zx=dudx, zy=dudy)
end

"""
    ZernikeBW(elements, ap, mask, norms)

Contains Zernike basis (in Born&Wolf norming) with the aperture, plotting mask, and element norms. Can be constructed as
    `ZernikeBW(dom::CartesianDomain2D, d::Real, maxorder::Integer)`,
    where `d` is the aperture *diameter* (not radius)
and
    `ZernikeBW(gridsize::Integer, maxorder::Integer)`

"""
struct ZernikeBW <: OrthogonalBasis
    elements::VectorOfArray
    ap::Array
    mask::Array
    norms::Vector
    function ZernikeBW(elements, ap, mask)
        return new(elements, ap, mask, [sqrt.(inner(f, ap .* f)) for f in elements])
    end
end

function ZernikeBW(gridsize::Integer, maxorder::Integer)
    return ZernikeBW(makezerniketable(gridsize, maxorder), makeaperture(gridsize)...)
end
function ZernikeBW(dom::CartesianDomain2D, d::Real, maxorder::Integer)
    return ZernikeBW(makezerniketable(dom, maxorder, d / 2), aperture(dom, d)...)
end

# function makezerniketable(gridsize::Integer, maxorder::Integer)
#     x = range(-1, 1, length=gridsize)
#     y = range(-1, 1, length=gridsize)
#     totalznum = Int((maxorder + 2) * (maxorder + 1) / 2)
#     ztable = [zernike(xc, yc, maxorder)[:z] for xc ∈ x,  yc ∈ y ]
#     zvec = VectorOfArray([zeros(gridsize, gridsize) for i = 1:totalznum])
#     [zvec[i,:] = ztable[i] for i = eachindex(ztable) ]
#     return zvec
# end

function makezerniketable(dom::CartesianDomain2D, maxorder::Integer, scale=1)
    x = dom.xrange / scale
    y = dom.yrange / scale
    totalznum = Int((maxorder + 2) * (maxorder + 1) / 2)
    ztable = [zernike(xc, yc, maxorder)[:z] for yc in y, xc in x]
    zvec = VectorOfArray([zeros(length(y), length(x)) for i in 1:totalznum])
    [zvec[i, :] = ztable[i] for i in eachindex(ztable)]
    return zvec
end

function makezerniketable(gridsize::Integer, maxorder::Integer)
    return makezerniketable(
        CartesianDomain2D(range(-1, 1; length=gridsize), range(-1, 1; length=gridsize)),
        maxorder,
    )
end

function makeaperture(gridsize::Integer, δ=0.0)
    x = range(-1, 1; length=gridsize)
    y = range(-1, 1; length=gridsize)
    # δ = 0. # tuning of the aperture size
    r = 1 + δ / gridsize
    ap = [(xc^2 + yc^2) <= r^2 ? 1 : 0 for xc in x, yc in y]
    # area = +(ap[:]...)
    phmask = [(xc^2 + yc^2) <= r^2 ? 1 : NaN for xc in x, yc in y]
    return (ap, phmask)
end

## Numbering schemes

triangle(k::Int) = k * (k + 1) ÷ 2
maxtriangle(i) = floor(Int, sqrt(2i + 0.25) - 0.5)

"""
    osa_j_to_nm(j::Int)

Convert OSA/ANSI standard indexing to (n,m) pair. Return named tuple.

```jldoctest
julia> osa_j_to_nm(0)
(n = 0, m = 0)

julia> osa_j_to_nm(65)
(n = 10, m = 10)

julia> osa_j_to_nm(495)
(n = 30, m = 30)
```
"""
function osa_j_to_nm(j::Int)
    n = maxtriangle(j)
    m = 2(j - triangle(n)) - n
    return (n=n, m=m)
end

"""
    nm_to_osa_j(;n,m)

Convert (n,m) index of Zernike polynomial to OSA/ANSI standard indexing.

```jldoctest
julia> nm_to_osa_j(n=4,m=0)
12

julia> nm_to_osa_j(m=4,n=4)
14

julia> nm_to_osa_j(n=30, m=30)
495

```
"""
function nm_to_osa_j(; n::Int, m::Int)
    j = (m + n) ÷ 2 + triangle(n)
    return j
end

export osa_j_to_nm, nm_to_osa_j

## make basis callable for convenience
# This is pixelated, the modal phase is better
# (basis::ZernikeBW)(j::Int) = elements(basis)[j + 1]
# (basis::ZernikeBW)(; n::Int, m::Int) = (basis)(nm_to_osa_j(; n=n, m=m))
(basis::ZernikeBW)(j::Int) = ModalPhase([j + 1], [1.0], basis)
function (basis::ZernikeBW)(; n::Int, m::Int)
    return ModalPhase([nm_to_osa_j(; n=n, m=m) + 1], [1.0], basis)
end

# Sparse Zernike
struct ZernikeBWSparse <: OrthogonalBasis
    elements::VectorOfArray
    ap::SparseMatrixCSC{Float64,Int64}
    mask::SparseMatrixCSC{Bool,Int64}
    norms::Vector
end

function makesparsezerniketable(dom::CartesianDomain2D, maxorder::Integer, apD=2, scale=1)
    x = dom.xrange / scale
    y = dom.yrange / scale
    ap = sparse(@. x^2 + y'^2 <= apD^2 / 4)
    is, js, vs = findnz(ap)
    totalznum = Int((maxorder + 2) * (maxorder + 1) / 2)
    zvec = VectorOfArray([similar(ap, Float64) for k in 1:totalznum])
    for ind in eachindex(is)
        zvec[is[ind], js[ind], :] .= zernike(x[is[ind]], y[js[ind]], maxorder)[:z]
    end
    norms = [norm(zer.nzval) for zer in zvec]
    norms /= norms[1]
    return zvec, ap, ap, norms
end

function makesparsezerniketable(gridsize::Integer, maxorder::Integer, apD=2)
    return makesparsezerniketable(
        CartesianDomain2D(range(-1, 1; length=gridsize), range(-1, 1; length=gridsize)),
        maxorder,
        apD,
    )
end

function ZernikeBWSparse(gridsize::Integer, maxorder::Integer)
    zer, ap, mask, norms = makesparsezerniketable(gridsize, maxorder)

    return ZernikeBWSparse(zer, ap, mask, norms)
end
function ZernikeBWSparse(dom::CartesianDomain2D, d::Real, maxorder::Integer)
    zer, ap, mask, norms = makesparsezerniketable(dom, maxorder, d)
    return ZernikeBWSparse(zer, ap, mask, norms)
end
