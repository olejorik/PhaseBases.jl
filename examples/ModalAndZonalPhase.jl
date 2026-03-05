using PhaseBases
import PhaseBases: decompose, decompose_and_complement, project
using LinearAlgebra: norm
using CairoMakie
CairoMakie.activate!(; type="png")

# # Modal and Zonal Phases
#
# PhaseBases represents wavefronts in two complementary forms:
#
# - **`ModalPhase`** — coefficient vector tied to a basis (lazy, no grid until `collect`)
# - **`ZonalPhase`** — explicit 2D pixel array
#
# This tutorial covers construction, arithmetic, and conversion between them.

# ## 1 — ModalPhase Construction

zbas = ZernikeBW(128, 6)

# ### Dense: full coefficient vector

coef = zeros(length(zbas))
coef[5] = 0.8; coef[13] = -0.4
wf = ModalPhase(coef, zbas)

# ### Sparse: only the non-zero terms

wf_sp = ModalPhase([5, 13], [0.8, -0.4], zbas)  ## same result

collect(wf) ≈ collect(wf_sp)   ## true

# ### Zero-initialised

wf0 = ModalPhase(zbas)         ## all-zero coefficients

# Access and mutate coefficients:

coefficients(wf)[5]            ## 0.8
coefficients!(wf0, coef)       ## set in place

# ## 2 — Arithmetic
#
# `ModalPhase` supports `+`, `-`, scalar `*` when phases share the same basis.

wf_a = ModalPhase([5],  [1.0], zbas)
wf_b = ModalPhase([13], [0.5], zbas)

wf_sum = wf_a + wf_b
coefficients(wf_sum)[[5, 13]]   ## [1.0, 0.5]

wf_neg = -wf_a
wf_scl = 3.0 * wf_b

# `norm` returns the coefficient-weighted norm (valid for orthogonal bases):

norm(wf_a)

# ## 3 — ZonalPhase
#
# A thin wrapper around a 2D `Float64` array.

arr = collect(wf)
zp = ZonalPhase(arr)

coefficients(zp) === arr        ## same underlying data
zp[64, 64]                      ## pixel access

# Arithmetic works element-wise:

zp2 = zp + zp
maximum(abs, coefficients(zp2) .- 2 .* arr)  ## ≈ 0

# ## 4 — Conversions
#
# **Modal → Zonal**: materialise via `ZonalPhase(modal)` or `convert`.

zp_from_modal = ZonalPhase(wf)

# **Zonal → Modal**: requires a decomposition step.

fitted = decompose(coefficients(zp), zbas)
wf_roundtrip = ModalPhase(fitted, zbas)

fig = Figure(; size=(500, 220))
ax1 = Axis(fig[1, 1]; title="original", aspect=DataAspect())
ax2 = Axis(fig[1, 2]; title="round-tripped", aspect=DataAspect())
heatmap!(ax1, collect(wf) .* mask(zbas); colormap=:RdBu)
heatmap!(ax2, collect(wf_roundtrip) .* mask(zbas); colormap=:RdBu)
fig

# ## 5 — Decompose and Complement
#
# `decompose_and_complement` splits an array into the in-basis projection
# and the orthogonal residual — useful for checking fit quality.

## add some high-frequency noise that the basis cannot represent
noisy = collect(wf) .+ 0.05 .* randn(128, 128) .* aperture(zbas)

fitted_coef, residual = decompose_and_complement(noisy, zbas)

fig2 = Figure(; size=(650, 220))
ax1 = Axis(fig2[1, 1]; title="input", aspect=DataAspect())
ax2 = Axis(fig2[1, 2]; title="projection", aspect=DataAspect())
ax3 = Axis(fig2[1, 3]; title="residual", aspect=DataAspect())
heatmap!(ax1, noisy .* mask(zbas); colormap=:RdBu)
heatmap!(ax2, project(noisy, zbas) .* mask(zbas); colormap=:RdBu)
heatmap!(ax3, residual .* mask(zbas); colormap=:RdBu)
fig2

# ## Summary
#
# | Feature / function | Purpose |
# |:---|:---|
# | `ModalPhase(coef, basis)` | Lazy coefficient-based phase |
# | `ModalPhase(inds, coef, basis)` | Sparse constructor |
# | `coefficients(ph)` / `coefficients!(ph, c)` | Get / set coefficients |
# | `collect(ph)` | Evaluate on grid → 2D array |
# | `norm(ph)` | Coefficient-weighted norm |
# | `ZonalPhase(array)` | Pixel-value phase wrapper |
# | `ZonalPhase(modal)` | Convert modal → zonal |
# | `decompose(arr, basis)` | Fit array → coefficient vector |
# | `decompose_and_complement(arr, basis)` | Coefficients + residual |
# | `project(arr, basis)` | Compose the fitted coefficients |
