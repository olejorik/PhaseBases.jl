using DrWatson
@quickactivate "PhaseBases"
using PhaseBases
using Test
using LinearAlgebra
using RecursiveArrayTools

@testset "Zernike generation" begin
    @test zernike(0.5, 0.5, 3) == (
        z=[1.0, 0.5, 0.5, 0.5, 0.0, 0.0, 0.25, -0.25, -0.25, -0.25],
        zx=[0.0, 0.0, 1.0, 1.0, 2.0, 1.0, 1.5, 1.5, 1.0, 0.0],
        zy=[0.0, 1.0, 0.0, 1.0, 2.0, -1.0, 0.0, 1.0, 1.5, -1.5],
    )
    z = ZernikeBW(5, 4)
    @test elements(z, 15) == [
        -4.0 -0.4375 1.0 -0.4375 -4.0
        -0.4375 -0.25 0.0625 -0.25 -0.4375
        1.0 0.0625 0.0 0.0625 1.0
        -0.4375 -0.25 0.0625 -0.25 -0.4375
        -4.0 -0.4375 1.0 -0.4375 -4.0
    ]
    @test elements(z, 5) == [
        3.0 1.5 1.0 1.5 3.0
        1.5 0.0 -0.5 0.0 1.5
        1.0 -0.5 -1.0 -0.5 1.0
        1.5 0.0 -0.5 0.0 1.5
        3.0 1.5 1.0 1.5 3.0
    ]
    @test aperturedelements(z, 1) == [
        0.0 0.0 1.0 0.0 0.0
        0.0 1.0 1.0 1.0 0.0
        1.0 1.0 1.0 1.0 1.0
        0.0 1.0 1.0 1.0 0.0
        0.0 0.0 1.0 0.0 0.0
    ]

    coef = zeros(length(z))
    coef[10:11] = [-1, 1]
    ph = compose(z, coef)
    @test ph == [
        -2.0 -2.875 0.0 2.875 2.0
        1.75 -0.25 0.0 0.25 -1.75
        1.0 0.125 0.0 -0.125 -1.0
        -1.25 -0.25 0.0 0.25 1.25
        -2.0 0.125 0.0 -0.125 2.0
    ]
end

@testset "inner product" begin
    z = ZernikeBW(5, 4)
    el = z.elements
    @test PhaseBases.inner(el[1], el[1]) == 25.0
    @test PhaseBases.inner(el[1], el[2]) == 0

    a = reshape(1:27, (3, 3, 3))
    b = VectorOfArray([a[:, :, i] for i in 1:last(size(a))])
    coef = [1, 100, 10000]
    @test PhaseBases.inner(b, coef) ==
        PhaseBases.inner(coef, b) ==
        [
            191001 221304 251607
            201102 231405 261708
            211203 241506 271809
        ]
end

@testset "aperture" begin
    dom1 = PhaseBases.CartesianDomain2D(-1:0.2:1, -0.9:0.2:0.9)
    ap, mask = PhaseBases.aperture(dom1, 1)
    ap, mask = PhaseBases.aperture(dom1, 1)
    @test ap == [
        0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 1 0 0 0 0 0
        0 0 0 1 1 1 1 1 0 0 0
        0 0 0 1 1 1 1 1 0 0 0
        0 0 0 1 1 1 1 1 0 0 0
        0 0 0 1 1 1 1 1 0 0 0
        0 0 0 0 0 1 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0
    ]
end


@testset "decomposition" begin
    # The sampled basis is not orthogonal by default!
    zbig = ZernikeBW(512, 10)
    # m = PhaseBases.innermatrix(zbig);
    # norms = sqrt.(diag(m));
    # [zbig.elements[i] *= 1 / norms[i] for i in eachindex(zbig.elements)];
    # m = PhaseBases.innermatrix(zbig);
    # using ZChop
    # m = zchop(m);
    # r = m - I(length(zbig.elements));
    # maximum(abs.(r))

    length(zbig) == 66
    coefs = rand(length(zbig))
    aberration = compose(zbig, coefs)
    rest_coefs = decompose(aberration, zbig)
    @test all(rest_coefs .≈ coefs)
end

@testset "Zernike composition-decomposition consistency" begin
    ## test zero coefficients produce zero phase
    z = ZernikeBW(10, 3)
    coef_zero = zeros(length(z))
    ph_zero = compose(z, coef_zero)
    @test all(ph_zero .== 0.0)

    ## test single-mode composition and decomposition
    coef_single = zeros(length(z))
    idx = 4
    coef_single[idx] = 3.14
    ph_single = compose(z, coef_single)
    @test ph_single ≈ coef_single[idx] * elements(z, idx)
    rec_single = decompose(ph_single, z)
    @test rec_single ≈ coef_single

    ## test random coefficients round-trip
    coef_rand = randn(length(z))
    ph_rand = compose(z, coef_rand)
    rec_rand = decompose(ph_rand, z)
    @test all(rec_rand .≈ coef_rand)
end

@testset "API convenience methods" begin
    z = ZernikeBW(6, 3)

    ## test getindex shortcuts
    @test z[2] == elements(z, 2)
    sub = z[1:3]
    @test length(sub) == 3 && sub[1] == elements(z, 1)

    ## test dualelements API
    d = dualelements(z)
    @test length(d) == length(z)
    @test d[1] == getfield(z, :dualelements)[1]

    ## test indexes API
    idxs = indexes(z)
    @test isa(idxs, AbstractVector)
    @test all(x -> isa(x, Tuple) && length(x) == 2, idxs)
    @test idxs == getfield(z, :indexes)
end
