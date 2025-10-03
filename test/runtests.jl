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
    @test all(x -> isa(x, CartesianIndex) && length(x) == 2, idxs)
    @test idxs == getfield(z, :indexes)
end

@testset "aperturedelements mask option" begin
    z = ZernikeBW(5, 4)
    ## mask interface returns struct field (uses `isequal` for NaN checks)
    @test isequal(mask(z), getfield(z, :mask))

    ## default behavior (use_mask=false) matches original aperturedelements
    a1 = aperturedelements(z, 1)
    a_def = aperturedelements(z, 1; use_mask=false)
    @test isequal(a1, a_def)

    ## use_mask=true uses mask(z)
    m = mask(z)
    e1 = elements(z, 1)
    @test isequal(aperturedelements(z, 1; use_mask=true), m .* e1)

    ## multiple indices with use_mask
    inds = [1, 2]
    out = aperturedelements(z, inds; use_mask=true)
    @test length(out) == length(inds)
    @test all(isequal(out[i], m .* elements(z, inds[i])) for i in eachindex(inds))
end

@testset "PixelBasis construction" begin
    # Test construction from boolean aperture
    aperture = [
        false false false false false
        false true true true false
        false true true true false
        false true true true false
        false false false false false
    ]
    pb = PixelBasis(aperture)
    @test pb isa PixelBasis{Bool,2}
    @test size(pb.ap) == (5, 5)
    @test length(pb) == 9  # 3x3 inner square
    @test pb.size == (5, 5)

    # Test construction from Float64 aperture
    float_aperture = Float64.(aperture)
    pb_float = PixelBasis(float_aperture)
    @test pb_float isa PixelBasis{Float64,2}
    @test length(pb_float) == 9

    # Test construction from BitArray
    bit_aperture = BitArray(aperture)
    pb_bit = PixelBasis(bit_aperture)
    @test pb_bit isa PixelBasis{Bool,2}
    @test length(pb_bit) == 9

    # Test construction from explicit indexes
    manual_indexes = [CartesianIndex(2, 2), CartesianIndex(2, 3), CartesianIndex(3, 3)]
    pb_manual = PixelBasis(manual_indexes, (4, 4))
    @test length(pb_manual) == 3
    @test pb_manual.size == (4, 4)
end

@testset "PixelBasis elements interface" begin
    # Create simple 3x3 aperture for testing
    aperture = [
        true true true
        true false true
        true true true
    ]
    pb = PixelBasis(aperture)
    @test length(pb) == 8  # All pixels except center

    # Test single element access
    elem1 = elements(pb, 1)
    @test size(elem1) == (3, 3)
    @test sum(elem1) == 1  # Single unit element
    @test sum(.!iszero.(elem1)) == 1  # Only one non-zero element

    # Test multiple element access
    elems_multi = elements(pb, [1, 3, 5])
    @test length(elems_multi) == 3
    @test all(sum(e) == 1 for e in elems_multi)

    # Test bounds checking
    @test_throws BoundsError elements(pb, 0)
    @test_throws BoundsError elements(pb, length(pb) + 1)

    # Test lazy elements generator
    all_elems = elements(pb)
    count = 0
    for elem in all_elems
        count += 1
        @test sum(elem) == 1
    end
    @test count == length(pb)

    # Test collect all elements
    collected = collect(elements(pb))
    @test length(collected) == length(pb)
    @test all(sum(e) == 1 for e in collected)
end

@testset "PixelBasis compose and decompose" begin
    # Create test aperture
    aperture = [
        false true false
        true true true
        false true false
    ]
    pb = PixelBasis(aperture)
    @test length(pb) == 5  # Cross pattern

    # Test compose
    test_coef = [1.0, 2.0, 3.0, 4.0, 5.0]
    composed = compose(pb, test_coef)
    @test size(composed) == (3, 3)
    @test sum(.!iszero.(composed)) == 5
    @test sum(composed) == sum(test_coef)

    # Test decompose
    decomposed = decompose(composed, pb)
    @test length(decomposed) == 5
    @test decomposed ≈ test_coef

    # Test round-trip accuracy
    random_coef = randn(length(pb))
    roundtrip = decompose(compose(pb, random_coef), pb)
    @test roundtrip ≈ random_coef atol = 1e-14

    # Test error handling
    @test_throws ArgumentError compose(pb, [1.0, 2.0])  # Wrong length
    @test_throws ArgumentError decompose(zeros(2, 2), pb)  # Wrong size
end

@testset "PixelBasis orthonormality" begin
    # Create larger aperture for testing
    aperture = zeros(Bool, 5, 5)
    aperture[2:4, 2:4] .= true  # 3x3 block
    pb = PixelBasis(aperture)

    # Test norms are all 1 (orthonormal)
    norms_vec = norms(pb)
    @test length(norms_vec) == length(pb)
    @test all(norms_vec .≈ 1.0)

    # Test orthogonality of basis functions
    elem1 = elements(pb, 1)
    elem2 = elements(pb, 2)
    @test sum(elem1 .* elem2) == 0  # Orthogonal (disjoint support)

    # Test different elements have disjoint support
    all_elems = collect(elements(pb))
    for i in 1:length(pb), j in (i+1):length(pb)
        @test sum(all_elems[i] .* all_elems[j]) == 0
    end
end

@testset "PixelBasis interface compliance" begin
    # Create test basis
    test_aperture = [
        true false true
        false true false
        true false true
    ]
    pb = PixelBasis(test_aperture)

    # Test AbstractBasis interface methods
    @test aperture(pb) == test_aperture
    @test length(pb) == 5
    @test indexes(pb) isa Vector{CartesianIndex{2}}
    @test length(indexes(pb)) == 5

    # Test string representation
    str_repr = string(pb)
    @test contains(str_repr, "PixelBasis")
    @test contains(str_repr, "5 pixels")
end

@testset "PixelBasis with ModalPhase" begin
    # Create test aperture
    aperture = [
        true true
        true true
    ]
    pb = PixelBasis(aperture)

    # Test ModalPhase creation
    coef = [1.0, 2.0, 3.0, 4.0]
    mp = ModalPhase(coef, pb)
    @test mp isa ModalPhase
    @test length(mp.coef) == 4
    @test mp.basis === pb

    # Test collect (composition)
    phase_array = collect(mp)
    @test size(phase_array) == (2, 2)
    @test phase_array ≈ compose(pb, coef)

    # Test collect! (in-place composition)
    target = zeros(2, 2)
    result = collect!(target, mp)
    @test result === target  # Returns same array
    @test target ≈ phase_array  # Same result as collect
    @test target ≈ compose(pb, coef)

    # Test collect! with pre-existing values (should overwrite)
    target .= 999.0
    collect!(target, mp)
    @test target ≈ phase_array

    # Test collect! with different coefficient values
    new_coef = [10.0, 20.0, 30.0, 40.0]
    coefficients!(mp, new_coef)
    collect!(target, mp)
    @test target ≈ compose(pb, new_coef)

    # Test round-trip through ModalPhase
    reconstructed = decompose(collect(mp), pb)
    @test reconstructed ≈ new_coef atol = 1e-14

    # Test coefficient access and modification
    @test coefficients(mp) == new_coef
end

@testset "collect! with Zernike ModalPhase" begin
    # Test collect! with Zernike basis
    z = ZernikeBW(8, 3)
    coef = randn(length(z))
    mp_zern = ModalPhase(coef, z)

    # Test collect! allocation-free composition
    target = zeros(size(aperture(z)))
    result = collect!(target, mp_zern)
    @test result === target
    @test target ≈ collect(mp_zern)
    @test target ≈ compose(z, coef)

    # Test collect! with updated coefficients (auto-updating)
    coef .*= 2.0  # Modify original coefficient vector
    collect!(target, mp_zern)
    @test target ≈ compose(z, coef)

    # Test size compatibility check
    wrong_size_target = zeros(5, 5)
    @test_throws ArgumentError collect!(wrong_size_target, mp_zern)
end

@testset "PixelBasis performance characteristics" begin
    # Create larger aperture for performance testing
    aperture = zeros(Bool, 20, 20)
    center = 10.5
    for i in 1:20, j in 1:20
        if sqrt((i - center)^2 + (j - center)^2) <= 8
            aperture[i, j] = true
        end
    end
    pb = PixelBasis(aperture)
    n_pixels = length(pb)
    @test n_pixels > 100  # Reasonably large

    # Test compose performance (should be O(n_pixels))
    coef = randn(n_pixels)
    composed = compose(pb, coef)
    @test sum(.!iszero.(composed)) == n_pixels

    # Test decompose performance (should be O(n_pixels)) 
    decomposed = decompose(composed, pb)
    @test length(decomposed) == n_pixels
    @test decomposed ≈ coef atol = 1e-14

    # Test that lazy elements don't consume excessive memory
    elements_gen = elements(pb)
    @test elements_gen isa Base.Generator

    # Test iterator correctness without storing all
    element_count = 0
    sum_of_sums = 0.0
    for elem in elements_gen
        element_count += 1
        sum_of_sums += sum(elem)
    end
    @test element_count == n_pixels
    @test sum_of_sums ≈ n_pixels  # Each element sums to 1
end
