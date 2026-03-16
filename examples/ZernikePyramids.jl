# ```@meta
# CurrentModule = PhaseBases
# DocTestSetup = quote
#     using PhaseBases
# end
# ```

using PhaseBases
using CairoMakie
CairoMakie.activate!(; type="png")

# # Zernike Pyramids: Visualizing Ordering Conventions
#
# Zernike polynomials are traditionally displayed by arranging them in a table
# indexed by radial order $n$ (rows) and azimuthal frequency $m$ (columns).
# In the Born & Wolf / OSA convention this gives a triangular "pyramid."
#
# However, the **Fringe** (University of Arizona) convention re-orders the
# polynomials so that terms with the same optical path difference are grouped
# together. When arranged by their $(n, m)$ pairs, the Fringe ordering traces
# a **diamond**-shaped pa
# ttern rather than a triangle.
#
# This page renders both layouts side by side.

# ## Setup
#
# We create a moderately sized Zernike basis and a helper function that
# plots a single polynomial into a subplot cell.

zbas = ZernikeBW(128, 10)  ## 128×128 grid, polynomials up to radial order 10
ap = mask(zbas)            ## NaN outside aperture — clean display

"""
    plot_zernike!(fig, row, col, osa_idx; kw...)

Place a `heatmap` of the Zernike polynomial at 1-based array position
`osa_idx` (= OSA index + 1) into `fig[row, col]`.
"""
function plot_zernike!(fig, row, col, osa_idx; label="", kw...)
    ax = Axis(
        fig[row, col]; aspect=DataAspect(), title=label, titlesize=10, titlefont=:regular
    )
    heatmap!(
        ax,
        (elements(zbas, osa_idx) .* ap)';  ## transpose to match display orientation
        colormap=reverse(cgrad(:RdBu)),
        colorrange=(-1, 1),
        kw...,
    )
    hidedecorations!(ax)
    hidespines!(ax)
    return ax
end

# ## Born & Wolf / OSA Pyramid
#
# In the OSA convention, polynomials are ordered by increasing $n$, and within
# each order by increasing $m$ (from $-n$ to $+n$ in steps of 2). Laying them
# out on a grid of $(n, m)$ gives the classic triangle:
#
# ```
#             (0,0)
#           (1,-1)(1,1)
#         (2,-2)(2,0)(2,2)
#       (3,-3)(3,-1)(3,1)(3,3)
#     …
# ```

maxn = 6   ## show orders 0 through 6

fig_bw = Figure(; size=(900, 900))

for n in 0:maxn
    for m in (-n):2:n
        osa_j = nm_to_osa_j(; n=n, m=m)
        arr_idx = osa_j + 1   ## 1-based position in ZernikeBW
        ## Column: map m ∈ [-maxn, maxn] to grid columns 1..(2maxn+1)
        gcol = m + maxn + 1
        ## Row: n+1 (offset by 1 for the title row)
        grow = n + 2
        plot_zernike!(fig_bw, grow, gcol, arr_idx; label="($n,$m)")
    end
end

## Force equal cell sizes so isolated extremal columns (e.g. (6,-6)) are not oversized
cell_px = 65
for i in 1:(2maxn + 1)
    colsize!(fig_bw.layout, i, Fixed(cell_px))
end
for i in 2:(maxn + 2)
    rowsize!(fig_bw.layout, i, Fixed(cell_px))
end
Label(fig_bw[1, :], "Born & Wolf / OSA Pyramid  (n = 0 … $maxn)"; fontsize=18, font=:bold)
Colorbar(
    fig_bw[maxn + 3, :];
    colormap=reverse(cgrad(:RdBu)),
    limits=(-1, 1),
    vertical=false,
    width=Relative(0.5),
    height=Relative(0.05),
    label="Zernike Value",
)
resize_to_layout!(fig_bw)

fig_bw

# Each row has $n+1$ polynomials, forming the familiar triangular pyramid.
# The top is piston $(0,0)$, and each subsequent row adds one more
# azimuthal frequency on each side.


# ## Fringe / University of Arizona Diamond
#
# The Fringe ordering groups polynomials by "rings" of roughly equal optical
# significance. When we plot them by their $(n, m)$ coordinates, they trace a
# diamond (or rotated square) rather than a triangle, because the Fringe
# convention interleaves different radial orders.
#
# Below we show the first 36 Fringe terms (indices 1–36). With a basis of order 10,
# all 36 terms — including the bottom-apex terms F34=(9,-1), F35=(9,1), F36=(10,0) —
# are fully rendered.

n_fringe = 36  ## number of Fringe terms to display

## Collect all (n,m) pairs for Fringe 1..n_fringe
fringe_nm = [fringe_j_to_nm(j) for j in 1:n_fringe]
max_n_fringe = maximum(nm.n for nm in fringe_nm)
max_m_fringe = maximum(abs(nm.m) for nm in fringe_nm)

fig_fr = Figure(; size=(900, 1000))


for j in 1:n_fringe
    nm = fringe_nm[j]
    n, m = nm.n, nm.m
    osa_j = nm_to_osa_j(; n=n, m=m)
    arr_idx = osa_j + 1
    ## Need a larger basis if max order exceeds what we built
    arr_idx > length(zbas) && continue
    ## Grid position: row = n, col = m (centered)
    grow = n + 2
    gcol = m + max_m_fringe + 1
    plot_zernike!(fig_fr, grow, gcol, arr_idx; label="F$j ($n,$m)")
end

Label(fig_fr[1, :], "Fringe Diamond  (first $n_fringe terms)"; fontsize=18, font=:bold)
Colorbar(
    fig_fr[max_n_fringe + 3, :];
    colormap=reverse(cgrad(:RdBu)),
    limits=(-1, 1),
    vertical=false,
    width=Relative(0.5),
    height=Relative(0.05),
    label="Zernike Value",
)

## Force equal cell sizes so isolated extremal columns (e.g. (5,-5)) are not oversized
cell_px_fr = 65
for i in 1:(2max_m_fringe + 1)
    colsize!(fig_fr.layout, i, Fixed(cell_px_fr))
end
for i in 2:(max_n_fringe + 3)
    rowsize!(fig_fr.layout, i, Fixed(cell_px_fr))
end


resize_to_layout!(fig_fr)
fig_fr

# The diamond shape is clearly visible: Fringe index 1 sits at the top (piston),
# and the "rings" expand diagonally until they reach maximum width at the n=5 row
# (containing the (5,±5) trefoil terms), then contract back to a single tip at
# F36=(10,0) — the secondary spherical — at the bottom apex.
# Unlike the BW pyramid, high-$n$ low-$|m|$ terms (like primary spherical
# F9=(4,0)) appear early, while high-$|m|$ low-$n$ terms of the same order
# (like F17=(4,±4)) appear much later.


# ## Side-by-Side Index Maps
#
# As a compact reference, the following table shows how the first 28
# polynomials are numbered in each convention:

println(
    rpad("(n, m)", 10),
    " | ",
    rpad("OSA", 5),
    " | ",
    rpad("Noll", 5),
    " | ",
    rpad("Fringe", 6),
)
println("-"^35)
for n in 0:6
    for m in (-n):2:n
        osa = nm_to_osa_j(; n=n, m=m)
        noll = nm_to_noll_j(; n=n, m=m)
        fr = nm_to_fringe_j(; n=n, m=m)
        println(
            rpad("($n, $m)", 10),
            " | ",
            rpad(osa, 5),
            " | ",
            rpad(noll, 5),
            " | ",
            rpad(fr, 6),
        )
    end
end

# ## Summary
#
# - **OSA / Born & Wolf**: polynomials tile a **triangle** — each row $n$ has
#   $n+1$ terms, arranged symmetrically around $m=0$.
# - **Fringe**: polynomials tile a **diamond** — terms are grouped by optical
#   significance, interleaving different radial orders.
# - Both orderings describe the *same* polynomials; only the single-index
#   numbering differs. Use [`reorder`](@ref) or the `j_to_nm` / `nm_to_j`
#   dispatch to translate freely.
