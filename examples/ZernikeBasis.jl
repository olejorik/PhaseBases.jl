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
zbas(; n=3, m=1)           ## ModalPhase for coma

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

nm_to_fringe_j(; n=3, m=-1)   ## Fringe index for vertical coma
nm_to_osa_j(; n=3, m=-1)      ## OSA index for the same polynomial
nm_to_noll_j(; n=3, m=-1)     ## Noll index

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
coef0[5] = 1.0;
coef0[13] = 0.5;   ## tilt + spherical
arr = compose(zbas, coef0)

coef_dense = decompose(arr, zbas)
coef_sparse = decompose(arr, zbas_sp)

fig2 = Figure(; size=(500, 250))
ax = Axis(
    fig2[1, 1]; xlabel="basis index", ylabel="coefficient", title="decompose accuracy"
)
scatterlines!(ax, coef0; label="truth", markersize=8)
scatterlines!(ax, coef_dense; label="ZernikeBW", markersize=6)
scatterlines!(ax, coef_sparse; label="Sparse", markersize=6, linestyle=:dash)
axislegend(ax; position=:rt)
fig2

# ## 5 — Algorithm and Derivatives
#
# Both `ZernikeBW` and `ZernikeBWSparse` are backed by the efficient recurrence
# published by Andersen (2018):
#
# > T. B. Andersen, *"Efficient and robust recurrence relations for the Zernike
# > circle polynomials and their derivatives in Cartesian coordinates,"*
# > Opt. Express **26**(15), 18878, Jul. 2018.
# > <https://doi.org/10.1364/OE.26.018878>
#
# The key advantage of this recurrence is that it computes the exact Cartesian
# partial derivatives $\partial Z_j/\partial x$ and $\partial Z_j/\partial y$
# within the **same loop** as the polynomial values — no finite differences,
# no extra asymptotic cost.
#
# ### Low-level API
#
# The function `zernike(x, y, maxord)` evaluates all polynomials and their
# derivatives at a single point and returns a named tuple:

result = zernike(0.3, 0.4, 4)
## result.z  — polynomial values Z₁ … Z₁₅  (1-based, Andersen internal order)
## result.zx — ∂Z/∂x at (0.3, 0.4)
## result.zy — ∂Z/∂y at (0.3, 0.4)

# !!! note
#     `zernike` uses Andersen's internal ordering: within each radial order,
#     sine terms come first, index 1 = piston.  Use `nm_to_osa_j` etc. to
#     translate to your preferred convention.
#
# ### Derivative table
#
# The following figure shows $\partial Z_j / \partial y$ for all 66 polynomials
# up to radial order 10, rendered as filled contour plots — one panel per
# polynomial, laid out in a regular 6-column grid.  The structure mirrors the
# visualisations published at <https://olejorik.github.io/post/zernikecalculations/>.

const NCOLS_DER = 6
const MAXORD_DER = 10
const NZ_DER = Int((MAXORD_DER + 2) * (MAXORD_DER + 1) / 2)   ## 66
const NROWS_DER = cld(NZ_DER, NCOLS_DER)                       ## 11
const RES_DER = 128

## Evaluate ∂Z/∂y on a Cartesian grid; NaN outside the unit disk
xs_der = range(-1, 1; length=RES_DER)
zy_tab = [fill(NaN, RES_DER, RES_DER) for _ in 1:NZ_DER]
for (ci, xv) in enumerate(xs_der), (ri, yv) in enumerate(xs_der)
    xv^2 + yv^2 > 1 && continue
    vals = zernike(xv, yv, MAXORD_DER)
    for j in 1:NZ_DER
        zy_tab[j][ci, ri] = vals.zy[j]   ## [xi, yi] — Makie convention
    end
end

## Draw the grid: one contourf panel per polynomial
cell_px_der = RES_DER
fig_der = Figure(; size=(NCOLS_DER * cell_px_der + 20, NROWS_DER * cell_px_der + 20))
for j in 1:NZ_DER
    row = fld1(j, NCOLS_DER)
    col = mod1(j, NCOLS_DER)
    ax = Axis(fig_der[row, col]; aspect=DataAspect())
    contourf!(ax, xs_der, xs_der, zy_tab[j]; colormap=:Blues, levels=20)
    contour!(ax, xs_der, xs_der, zy_tab[j]; color=:white, levels=19)
    hidedecorations!(ax)
    hidespines!(ax)
end
resize_to_layout!(fig_der)
fig_der

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
# | `zernike(x, y, n)` | Low-level evaluation + exact Cartesian derivatives |
