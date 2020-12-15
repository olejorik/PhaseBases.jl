using DrWatson
@quickactivate "PhaseBases"
using PhaseBases
using Test
using LinearAlgebra
using RecursiveArrayTools

@testset "Zernike generation" begin
    @test zernike(.5, .5, 3) == (z = [1.0, 0.5, 0.5, 0.5, 0.0, 0.0, 0.25, -0.25, -0.25, -0.25], zx = [0.0, 0.0, 1.0, 1.0, 2.0, 1.0, 1.5, 1.5, 1.0, 0.0], zy = [0.0, 1.0, 0.0, 1.0, 2.0, -1.0, 0.0, 1.0, 1.5, -1.5])
    z = ZernikeBW(5, 4);
    @test z.elements[:,:,15] == 
   [-4.0     -0.4375  1.0     -0.4375  -4.0
    -0.4375  -0.25    0.0625  -0.25    -0.4375
     1.0      0.0625  0.0      0.0625   1.0
    -0.4375  -0.25    0.0625  -0.25    -0.4375
    -4.0     -0.4375  1.0     -0.4375  -4.0]
    @test PhaseBases.elements(z, [5,10])[1] == 
        [3.0 1.5 1.0 1.5 3.0; 1.5 0.0 -0.5 0.0 1.5; 1.0 -0.5 -1.0 -0.5 1.0; 1.5 0.0 -0.5 0.0 1.5; 3.0 1.5 1.0 1.5 3.0]
    @test  PhaseBases.aperturedelements(z, 1) == PhaseBases.aperturedelements(z, [1,2])[:,:,1] == [   
     0.0  0.0  1.0  0.0  0.0
     0.0  1.0  1.0  1.0  0.0
     1.0  1.0  1.0  1.0  1.0
     0.0  1.0  1.0  1.0  0.0
     0.0  0.0  1.0  0.0  0.0]


     coef = zeros(length(z))
     coef[10:11] = [-1,1]
     ph = compose(z, coef)
     @test ph == [
        -2.0   -2.875  0.0   2.875   2.0
        1.75  -0.25   0.0   0.25   -1.75
        1.0    0.125  0.0  -0.125  -1.0
       -1.25  -0.25   0.0   0.25    1.25
       -2.0    0.125  0.0  -0.125   2.0
     ]

end


@testset "inner product" begin
    z = ZernikeBW(5, 4);  
    el = z.elements;
    @test PhaseBases.inner(el[1],el[1]) == 25.0
    @test PhaseBases.inner(el[1], el[2]) == 0

    a = reshape(1:27,(3,3,3))
    b= VectorOfArray([a[:,:,i] for i =1:last( size(a))])
    @test PhaseBases.inner(b,[1,100,10000]) == PhaseBases.inner([1,100,10000],b) == [
        191001  221304  251607
        201102  231405  261708
        211203  241506  271809   ]

end

@testset "aperture" begin
dom1= PhaseBases.CartesianDomain2D(-1:.2:1, -.9:.2:.9)
ap, mask = PhaseBases.aperture(dom1, 1);
ap, mask = PhaseBases.aperture(dom1, 1);
@test ap == [
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
]
end

# The sampled basis is not orthogonal by default!
zbig = ZernikeBW(512, 10);
m = PhaseBases.innermatrix(zbig);
norms = sqrt.(diag(m));
[zbig.elements[i] *= 1 / norms[i] for i in eachindex(zbig.elements)];
m = PhaseBases.innermatrix(zbig);
using ZChop
m = zchop(m);
r = m - I(length(zbig.elements));
maximum(abs.(r))

