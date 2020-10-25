# PhaseBases

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://olejorik.github.io/PhaseBases.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://olejorik.github.io/PhaseBases.jl/dev)
[![Build Status](https://github.com/olejorik/PhaseBases.jl/workflows/CI/badge.svg)](https://github.com/olejorik/PhaseBases.jl/actions)
[![Coverage](https://codecov.io/gh/olejorik/PhaseBases.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/olejorik/PhaseBases.jl)

WIP

PhaseBases.jl is a small package implementing typical bases used for decomposition of the phase of the optical field in the pupil, like Zernike polynomials and 
Gaussian radial base functions.

The package is intended for the use scenario with iterative algorithms (like phase retrieval), where different combinations of the basis functions should be evaluated multiple times.
For this purpose, the basis function are precalculated on a fixed grid and kept for future access.