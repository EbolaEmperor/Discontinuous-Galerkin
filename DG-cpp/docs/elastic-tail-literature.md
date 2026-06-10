# Elastic Tail Literature Notes

Downloaded PDFs live in `docs/literature/`.

## Sources Used

- `2002_Peskin_immersed_boundary_method.pdf`: Peskin's Acta Numerica review of
  immersed elastic boundaries and Eulerian/Lagrangian coupling.
- `2011_Tian_Luo_Zhu_Liao_elastic_filaments_IB_LBM.pdf`: elastic filaments near
  a cylinder wake, including Karman-gait and entrainment configurations.
- `2003_Liao_Karman_gait_trout_cylinder_wake.pdf`: experimental motivation for
  station keeping in cylinder wakes.
- `2011_Heltai_Costanzo_variational_immersed_fem.pdf`: variational finite
  element immersed coupling and the adjoint relationship between interpolation
  and force spreading.

The Zhu-Peskin 2002 JCP flapping-filament article was read from the online PDF
view/search result, but the publisher mirror repeatedly failed via command-line
download. It is still cited in the implementation rationale because Tian et al.
explicitly reuse that filament force discretization.

## Implementation Consequences

- The elastic tail should be represented as a Lagrangian elastic boundary with
  mass and bending/stretching energy, not as a rendered polyline that samples
  local drag only.
- Fluid/structure transfer should use paired operators: interpolate velocity
  from the DG field to the tail and scatter the reaction force back into the
  DG RHS. The existing `IBCoupler` already implements this FE-IB structure.
- A moving head clamp must prescribe both position and velocity at the clamped
  tail nodes. Treating the clamp as stationary after teleporting it to the head
  injects nonphysical high-frequency energy.
- The tail should be almost inextensible.  After the implicit Cosserat step, the
  elastic tadpole projects segment lengths back to their rest values from the
  clamp toward the tip; this prevents axial stretch from dominating the intended
  bending response.
- Light immersed structures are sensitive to added-mass instability in explicit
  partitioned coupling. The new elastic tadpole therefore uses relaxed direct
  forcing, force capping, and a mild 1-2-1 smoothing pass along the tail.
- The original `navier_stokes_tadpole` should remain the rigid-tail drift demo.
  The elastic-tail model is a separate executable and config:
  `navier_stokes_tadpole_elastic` with `examples/ns_tadpole_elastic_config.json`.
