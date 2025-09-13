# Discontinuous-Galerkin
A Matlab implementation of Discontinuous Galerkin methods

- `Poisson1D` solves the following equation using $C^{-1}P_k\;(k\geq 1)$ elements.
  $$
  \left\{\begin{array}{l}
  -u''=f\\
  u(x_0)=g_0,\;u(x_1)=g_1
  \end{array}
  \right.
  $$

- `Poisson2D` solves the following equation using $C^{-1}P_k\;(k\geq 1)$ elements.
  $$
  \left\{\begin{array}{ll}
  -\Delta u=f,\quad &\text{in}\;\Omega\\
  u=g,\quad &\text{on}\;\partial \Omega
  \end{array}
  \right.
  $$

- 

- `Biharmonic1D` solves the following equation using $C^{0}P_k\;(k\geq 2)$ elements.
  $$
  \left\{\begin{array}{l}
  u^{(4)}=f\\
  u(x_0)=g_0,\;u(x_1)=g_1\\
  u'(x_0)=h_0,\;u'(x_1)=h_1
  \end{array}
  \right.
  $$

- `Biharmonic2D` solves the following equation using $C^{0}P_k\;(k\geq 1)$ elements.
  $$
  \left\{\begin{array}{ll}
  \Delta^2 u=f,\quad &\text{in}\;\Omega\\
  u=g,\quad &\text{on}\;\partial \Omega\\
  n\cdot \nabla u=h,\quad &\text{on}\;\partial \Omega
  \end{array}
  \right.
  $$

