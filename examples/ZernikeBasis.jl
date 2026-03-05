using PhaseBases
import PhaseBases: decompose
using CairoMakie
CairoMakie.activate!(; type="png")

# # Zernike Bases
#
# PhaseBases ships two Zernike basis types that store Born & Wolf–normalised
# polynomials on a computational grid.  Both share the same constructors
# but differ in storage and decomposition strategy.

# ## 1 — ZernikeBW: Dense Basis with Pseudo-Inverse
#
# `ZernikeBW` stores every polynomial as a full 2D array and computes
# a pseudo-inverse for accurate decomposition (no orthogonality assumption).

zbas = ZernikeBW(128, 6)
length(zbas)          ## 28 polynomials (orders 0–6)

# Two constructor forms:
# - `ZernikeBW(gridsize, maxorder)` — unit circle on a square grid
# - `ZernikeBW(dom, diameter, maxorder)` — arbitrary `CartesianDomain2D`

# ### Inspecting the structure

aperture(zbas) |> size   ## (128, 128) binary mask
norms(zbas)[1:5]         ## Born & Wolf norms (piston norm = √π ≈ 1.77)

# Elements are accessible by integer index or by (n, m):

zbas(4)                  ## ModalPhase for OSA j = 4 (defocus)
zbas(n=3, m=1)           ## ModalPhase for coma

# ### Visualising individual polynomials

fig = Figure(; size=(650, 280))
for (i, j) in enumerate([1, 5, 8, 13])
    ax = Axis(fig[1, i]; title="Z$j", aspect=DataAspect())
    heatmap!(ax, elements(zbas)[j] .* mask(zbas); colormap=:RdBu)
    hidedecorations!(ax)
end
fig

# ## 2 — ZernikeBWSparse: Sparse + Orthogonality Assumption
#
# `ZernikeBWSparse` stores polynomials as `SparseMatrixCSC` columns.
# It subtypes `OrthogonalBasis`, so `decompose` divides by norm²
# instead of using a pseudo-inverse — faster but approximate on
# pixelated grids.

zbas_sp = ZernikeBWSparse(128, 6)
typeof(elements(zbas_sp))   ## SparseMatrixCSC or VectorOfArray of sparse

# ### When to use which
#
# | | `ZernikeBW` | `ZernikeBWSparse` |
# |:---|:---|:---|
# | Storage | dense arrays | sparse matrices |
# | Decompose method | pseudo-inverse (accurate) | orthogonality (fast, approximate) |
# | Subtype | `AbstractBasis` | `OrthogonalBasis` |
# | Best for | moderate order, precise fits | high order, large grids |

# ## 3 — Index Conventions
#
# Six convention pairs convert between a single index `j` and the
# double index `(n, m)`. The standalone functions work on plain integers.

fringe_j_to_nm(4)    ## (n=2, m=0) — defocus in Fringe
noll_j_to_nm(4)      ## (n=2, m=0) — same polynomial, Noll j=4 too
osa_j_to_nm(4)       ## (n=2, m=2) — different polynomial in OSA!

# Inverse direction:

nm_to_fringe_j(n=3, m=-1)   ## Fringe index for vertical coma
nm_to_osa_j(n=3, m=-1)      ## OSA index for the same polynomial
nm_to_noll_j(n=3, m=-1)     ## Noll index

# Round-trip check:

all(nm_to_fringe_j(fringe_j_to_nm(j)) == j for j in 1:36)

# The `zerniketicks` helper generates tick labels for Makie plots:

ticks = zerniketicks(zbas)
ticks[1][1:5]   ## indices
ticks[2][1:5]   ## "(n, m)" labels

# ## 4 — Decompose / Compose Comparison
#
# Building a test wavefront and round-tripping through both basis types.

coef0 = zeros(length(zbas))
coef0[5] = 1.0; coef0[13] = 0.5   ## tilt + spherical
arr = compose(zbas, coef0)

coef_dense  = decompose(arr, zbas)
coef_sparse = decompose(arr, zbas_sp)

fig2 = Figure(; size=(500, 250))
ax = Axis(fig2[1, 1]; xlabel="basis index", ylabel="coefficient",
    title="decompose accuracy")
scatterlines!(ax, coef0; label="truth", markersize=8)
scatterlines!(ax, coef_dense; label="ZernikeBW", markersize=6)
scatterlines!(ax, coef_sparse; label="Sparse", markersize=6, linestyle=:dash)
axislegend(ax; position=:rt)
fig2

# ## Summary
#
# | Feature / function | Purpose |
# |:---|:---|
# | `ZernikeBW(N, order)` | Dense Zernike basis, accurate pseudo-inverse |
# | `ZernikeBWSparse(N, order)` | Sparse storage, fast orthogonal decompose |
# | `elements(b)`, `b[i]` | Access basis functions |
# | `norms(b)`, `aperture(b)`, `mask(b)` | Basis metadata |
# | `zbas(j)`, `zbas(n=, m=)` | Callable shortcut → `ModalPhase` |
# | `fringe_j_to_nm`, `noll_j_to_nm`, … | Index convention conversions |
# | `zerniketicks(b)` | Tick labels for plots |
# | `decompose(arr, b)` | Fit array to basis |
# | `compose(b, coef)` | Evaluate from coefficients |
