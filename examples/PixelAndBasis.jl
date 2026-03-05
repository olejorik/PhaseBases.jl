using PhaseBases
import PhaseBases: decompose, compose
using CairoMakie
CairoMakie.activate!(; type="png")

# # Pixel Basis and Generic Basis
#
# Besides Zernike polynomials, PhaseBases provides two general-purpose
# basis types for arbitrary function sets on a pixel aperture.

# ## 1 — PixelBasis
#
# `PixelBasis` is an orthonormal basis whose elements are unit-impulse
# functions at each aperture pixel. It is a thin, memory-efficient
# wrapper — elements are generated lazily on demand.

## simple circular aperture
N = 64
xs = range(-1, 1; length=N)
ap = [x^2 + y^2 ≤ 1.0 ? 1.0 : 0.0 for x in xs, y in xs]
pbas = PixelBasis(ap)
length(pbas)   ## number of non-zero pixels

# `compose` / `decompose` are O(N_pixels) direct-indexing operations:

## round-trip: array → coefficients → array
arr = randn(N, N) .* ap
coef = decompose(arr, pbas)
arr2 = compose(pbas, coef)
maximum(abs, (arr .- arr2) .* ap)   ## ≈ 0

# Alternative constructors:

pbas_bool = PixelBasis(ap .> 0 |> BitArray)   ## from BitArray
pbas_idx  = PixelBasis(indexes(pbas), (N, N)) ## from explicit indices

# ## 2 — Generic Basis (Arbitrary Functions)
#
# The `Basis` type builds a pseudo-inverse–based basis from *any* set
# of functions restricted to a known pixel subset.  This is the most
# general constructor in PhaseBases — use it when your functions are
# not Zernike polynomials.

## hand-craft three Gaussian bumps as basis functions
σ = 0.25
centers = [(-0.4, 0.0), (0.4, 0.0), (0.0, 0.5)]
funcs = [
    [exp(-((x - cx)^2 + (y - cy)^2) / (2σ^2)) for x in xs, y in xs]
    for (cx, cy) in centers
]

## restrict to the circular aperture
ap_idx = findall(!iszero, ap)
gbas = Basis(funcs, ap_idx)

length(gbas)         ## 3
norms(gbas)          ## Born & Wolf–style norms

# Visualise the three basis functions:

fig = Figure(; size=(550, 180))
for i in 1:3
    ax = Axis(fig[1, i]; title="g$i", aspect=DataAspect())
    heatmap!(ax, elements(gbas)[i] .* aperture(gbas); colormap=:viridis)
    hidedecorations!(ax)
end
fig

# ### Fitting a wavefront in the custom basis

## true wavefront: weighted sum of the Gaussians + noise
true_coef = [0.5, -0.3, 0.8]
wf_true = compose(gbas, true_coef) .+ 0.02 .* randn(N, N) .* ap

fitted = decompose(wf_true, gbas)
residual = wf_true .- compose(gbas, fitted)

fig2 = Figure(; size=(650, 220))
ax1 = Axis(fig2[1, 1]; title="input", aspect=DataAspect())
ax2 = Axis(fig2[1, 2]; title="fit", aspect=DataAspect())
ax3 = Axis(fig2[1, 3]; title="residual", aspect=DataAspect())
heatmap!(ax1, wf_true .* aperture(gbas); colormap=:RdBu)
heatmap!(ax2, compose(gbas, fitted) .* aperture(gbas); colormap=:RdBu)
heatmap!(ax3, residual .* aperture(gbas); colormap=:RdBu)
fig2

# Fitted coefficients vs truth:

hcat(true_coef, fitted)

# ## 3 — When to Use Which
#
# | Basis | Use case |
# |:---|:---|
# | `ZernikeBW` | Standard Zernike expansion, moderate order |
# | `ZernikeBWSparse` | High-order Zernike, large grids |
# | `PixelBasis` | Pixel-level operations, zonal phase |
# | `Basis` | Any custom function set (Gaussians, wavelets, …) |
# | `ShiftedBasis` | `Basis` with a non-zero origin (mean subtraction) |

# ## Summary
#
# | Feature / function | Purpose |
# |:---|:---|
# | `PixelBasis(aperture)` | Orthonormal pixel-impulse basis |
# | `PixelBasis(mask::BitArray)` | From boolean mask |
# | `Basis(functions, indexes)` | Generic pseudo-inverse basis |
# | `ShiftedBasis(funcs, origin, idx)` | Basis with offset origin |
# | `compose(b, coef)` | Coefficients → array |
# | `decompose(arr, b)` | Array → coefficients |
# | `elements(b)`, `norms(b)` | Inspect basis functions |
# | `mask(b)`, `aperture(b)` | Aperture metadata |
