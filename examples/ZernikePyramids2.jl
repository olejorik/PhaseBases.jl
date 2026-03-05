# ```@meta
# CurrentModule = PhaseBases
# DocTestSetup = quote
#     using PhaseBases
# end
# ```

using PhaseBases
using CairoMakie
CairoMakie.activate!(; type="png")

# # Zernike Pyramids: BBox Layout
#
# This page reproduces the BW pyramid and Fringe diamond from
# [ZernikePyramids](ZernikePyramids.md), but uses Makie's `BBox` positioning
# to place each subplot at an *exact pixel coordinate* derived
# analytically from the polynomial's $(n, m)$ indices.
#
# The mapping is simply:
# $$c_x = x_0 + m \cdot \tfrac{\delta}{2}, \quad c_y = y_{\text{top}} - n \cdot \delta$$
# where $\delta$ is the cell size in pixels and $(c_x, c_y)$ is the cell centre.
# Because every cell is placed independently, no empty grid cells exist and
# no `colsize!`/`rowsize!` corrections are needed.

# ## Setup

zbas = ZernikeBW(64, 10)  ## 64×64 grid, polynomials up to radial order 10
ap = mask(zbas)            ## NaN outside aperture

## Cell size, step, and padding (pixels).
## CELL  — width/height of each individual heatmap subplot.
## STEP  — centre-to-centre distance between adjacent cells (must be ≥ CELL).
##         gap between cells = STEP - CELL.
const CELL  = 64
const STEP  = 76   ## increase for more spacing, decrease toward CELL for tighter layout
const PAD_H = 50   ## horizontal padding on each side
const PAD_T = 60   ## space reserved for the title at the top
const PAD_B = 20   ## bottom padding

## Place a single Zernike heatmap centred at pixel (cx, cy) in fig.
function place_zernike!(fig, cx, cy, osa_arr_idx; label="")
    half = CELL / 2
    ax = Axis(fig;
        bbox       = BBox(cx - half, cx + half, cy - half, cy + half),
        aspect     = DataAspect(),
        title      = label,
        titlesize  = 9,
        titlefont  = :regular,
    )
    heatmap!(ax, elements(zbas, osa_arr_idx) .* ap; colormap=:RdBu)
    hidedecorations!(ax)
    hidespines!(ax)
    return ax
end


# ## Born & Wolf / OSA Pyramid
#
# Piston (0,0) sits at the apex; each subsequent row adds one term on each
# side, forming the classic triangle.
# With the BBox approach, the figure width and height are set
# analytically from `maxn`, so no post-layout resize step is needed.

maxn = 6   ## show radial orders 0 … 6

## Figure dimensions
bw_width  = (2maxn + 1) * STEP + 2PAD_H
bw_height = (maxn + 1)  * STEP + PAD_T + PAD_B

## Origin: horizontal centre; vertical top (below title padding)
bw_ox    = bw_width / 2
bw_top_y = bw_height - PAD_T - STEP / 2

fig_bw = Figure(; size=(bw_width, bw_height))

## Title using absolute scene coordinates
text!(fig_bw.scene,
    "Born & Wolf / OSA Pyramid  (n = 0 … $maxn)";
    position  = (bw_width / 2, bw_height - PAD_T / 2),
    align     = (:center, :center),
    fontsize  = 18,
    font      = :bold,
)

for n in 0:maxn
    for m in (-n):2:n
        arr_idx = nm_to_osa_j(n=n, m=m) + 1   ## 1-based position in zbas
        cx = bw_ox    + m * (STEP / 2)
        cy = bw_top_y - n * STEP
        place_zernike!(fig_bw, cx, cy, arr_idx; label="($n,$m)")
    end
end

fig_bw

# Because adjacent $m$ values within a row differ by 2,
# the horizontal step is $m \times \mathtt{STEP}/2$, placing them exactly $\mathtt{STEP}$
# apart. The gap between cell borders is $\mathtt{STEP} - \mathtt{CELL}$ pixels.


# ## Fringe / University of Arizona Diamond
#
# The same formula places Fringe terms on a diamond lattice: the layout
# reaches maximum width at $n = 5$ (the $\pm 5$ trefoil terms) and
# contracts back to a single tip at F36 = $(10, 0)$.

n_fringe = 36  ## number of Fringe terms to display

fringe_nm  = [fringe_j_to_nm(j) for j in 1:n_fringe]
max_n_fr   = maximum(nm.n for nm in fringe_nm)    ## 10
max_abs_m  = maximum(abs(nm.m) for nm in fringe_nm)  ## 5

## Figure dimensions
fr_width  = (2max_abs_m + 1) * STEP + 2PAD_H
fr_height = (max_n_fr + 1)   * STEP + PAD_T + PAD_B

fr_ox    = fr_width / 2
fr_top_y = fr_height - PAD_T - STEP / 2

fig_fr = Figure(; size=(fr_width, fr_height))

text!(fig_fr.scene,
    "Fringe Diamond  (first $n_fringe terms)";
    position  = (fr_width / 2, fr_height - PAD_T / 2),
    align     = (:center, :center),
    fontsize  = 18,
    font      = :bold,
)

for j in 1:n_fringe
    n, m   = fringe_nm[j].n, fringe_nm[j].m
    arr_idx = nm_to_osa_j(n=n, m=m) + 1
    arr_idx > length(zbas) && continue  ## guard (shouldn't trigger with order 10)
    cx = fr_ox    + m * (STEP / 2)
    cy = fr_top_y - n * STEP
    place_zernike!(fig_fr, cx, cy, arr_idx; label="F$j ($n,$m)")
end

fig_fr

# The diamond narrows smoothly because the Fringe convention includes
# only terms where the combination of $n$ and $|m|$ is below a threshold,
# naturally bounding the maximum azimuthal frequency at each radial order.
