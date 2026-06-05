# Discontinuous-Galerkin
A Matlab implementation of Discontinuous Galerkin methods

- `Poisson 1D/2D` solves the Poisson equation using $C^{-1}P_k\;(k\geq 1)$ elements.
- `Biharmonic 1D/2D` solves the Biharmonic equation using $C^{0}P_k\;(k\geq 2)$ elements.
- `PoissonMixed 1D/2D` solves the mixed-form Poisson equation using HDG method with $(C^{-1}P_k,C^{-1}P_k)$ mixed-element. And interpolate $u_h$ to $C^{-1}P_{k+1}$ space (by solving local Poisson equations) to show the super-convergence.

