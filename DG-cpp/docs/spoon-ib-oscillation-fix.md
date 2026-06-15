# Killing the immersed-boundary checkerboard in the stirred-soup (spoon) demo

## Symptom

The crescent-"spoon" stirring demo (`navier_stokes_spoon`) developed grid-scale
**checkerboard ("Gibbs") oscillations** in the vorticity: a stippled, element-to-element
sign-alternating texture generated at the moving blade, shed into the wake as discrete
packets, which then *persisted and locally amplified* (e.g. the dense stippled blob that
dominated the original frame ~180). The flow elsewhere was physically reasonable
(the paddle sheds a counter-rotating vortex pair), but the schlieren-like grid noise made
the frames look broken.

## Diagnosis

The IB was a **semi-implicit Schur-complement no-slip constraint**: Lagrangian markers
filling the *whole* blade footprint at spacing `dmark ≈ h`, with `u_h(X_k)=V_k` imposed
*pointwise* after the viscous solve (CG on `G = I H⁻¹ Iᵀ + εI`). Frame forensics
(high-pass + divergence diagnostics) plus a literature survey (9 papers deep-read; see
`docs/literature/ib-gibbs/SYNTHESIS.md`) attributed the noise to **five distinct
mechanisms**:

| ID | Mechanism | Key reference |
|----|-----------|---------------|
| M1 | IB force applied *after* the PPE → injects divergence near the body (the lagged `ppe_div_damping` is a band-aid; removing it blows up at t≈4) | Bao et al. DFIB 2017; Kallemov 2016; IBSE 2017 |
| M2 | Hard pointwise constraint = derivative kink → the P2 basis rings | IBSE 2017; Hester 2021 |
| M3 | Whole-footprint markers at `s≈h` → ill-conditioned Schur whose near-null modes are **high-frequency force modes** (multiplier jitter) | Kallemov 2016 |
| M4 | Body crosses element faces each step → on a *discontinuous* basis the pointwise interpolation jumps as a marker crosses a face → moving-body transfer jitter | Yang 2009; Bao-Kaye-Peskin 2016 |
| M5 | No slope limiter → shed grid-scale vorticity has **no fluid-side dissipation**, so packets persist/amplify | (gap; flagged by all) |

The divergence field made M1/M3 visible directly: the baseline showed a thick grid-scale
**divergence checkerboard ring** hugging the blade — the per-step source that the next
PPE has to chase.

The decisive experiment: coarsening the markers to `s≈2h` and switching to a smooth
volume kernel (M3+M4) — *with no other change* — **removed the shed wake checkerboard**
at t=13, while the vortex pair stayed intact, and it ran *faster* than the baseline
(IB-CG iterations dropped from ~7 to ~2; 29 markers vs 132).

## Fix (implemented)

All three remedies are off by default and composable; the spoon demo enables them via
`spoon_config.json`. They live in `NSIntegrator` (`src/navier_stokes/NavierStokes.{h,cpp}`).

### R2 — regularize the transfer  *(primary fix; addresses M3, M4, attenuates M2)*
- **Coarser markers**, `spoon_dmark ≈ 2h` (`0.032` for `h=0.016`): moves the Schur
  complement off the conditioning cliff (Kallemov). This alone is what removes the wake noise.
- **Mollified Wendland-C2 kernel transfer** (`ib_kernel_fac`, in units of `h`): the marker
  interpolation/spreading becomes a normalised volume integral of a C², non-negative radial
  kernel over the elements within radius `δ = ib_kernel_fac·h`, instead of a pointwise basis
  evaluation. Same kernel for read and spread (so `G` stays SPD); per-marker weights
  renormalised to sum to 1. Because the kernel spans several elements with smoothly-varying
  weights, the interpolation is **continuous in marker position** even as a marker crosses a
  (discontinuous) element face — directly killing the M4 jitter (Yang 2009's smoothing
  criterion, adapted to unstructured triangles).
  `setIBConstraint(..., kernelDelta)`; `IBConstraint::rows` holds the per-marker element rows.

### R1 — reprojection  *(addresses M1)*
- `ib_subiters` constraint↔projection sub-iterations: after the Schur correction, an extra
  PPE projection removes the divergence the IB force injected (reusing the pressure
  factorization), optionally re-applying the constraint. `ib_end_project` chooses whether the
  last sub-step is the projection (div-free) or the constraint (exact no-slip).

### ~~R3 — artificial viscosity~~ (REMOVED 2026-06-15)

**Ablation study conclusion: AV is iatrogenic.** A controlled experiment (4 runs: ctrl/no-AV,
hv_fac=0.002, 0.01, 0.05; all with R1+R2 active) showed that the **ctrl run (no AV, no HV, R1+R2
only) was completely clean** through the entire t=0–25 domain — including the blade-withdrawal
window (frame 240) where AV was originally motivated. Meanwhile the canonical config with
`av_beta=0.05` **produced** the very checkerboard it was supposed to suppress.

**Root cause**: the per-element blend `u ← (1-g)u + g·u_diff` introduces a *discontinuity at
element faces* (neighbouring elements get different `g` values), which creates a new inter-element
jump each step. This is a **positive feedback loop**: AV-induced face discontinuity → sensor fires
on the new jump → stronger blend → more face discontinuity. The diffusion-residual sensor was
supposed to be zero on smooth fields, but the per-element blend itself *creates* the non-smooth
content it detects.

All AV parameters (`av_beta`, `av_sensor_lo`, `av_sensor_hi`) have been removed from
`spoon_config.json`. The code in `NavierStokes.h` still supports them (returns immediately when
`avBeta==0`), so they can be re-enabled for other flows if needed — but for this demo they are
**harmful, not helpful**.

Dead-ends left off by default (with caveats in the code): a per-element **modal** filter
(`filter_strength`) only damps intra-element high modes and barely fires on this flow; global
**hyperviscosity** (`hv_fac`) also unnecessary — ctrl run proves R1+R2 alone suffice.

## Config knobs (spoon_config.json)

```
"spoon_dmark":   0.032,   // ~2h marker spacing (R2)
"ib_kernel_fac": 1.5,     // Wendland-C2 kernel radius in units of h (R2/M4)
"ib_subiters":   1,       // constraint<->projection reprojections (R1)
```

No AV/HV parameters needed. Reproduce with
`./build/navier_stokes_spoon spoon_config.json`.

## Result
Validated full t=42 run (441 frames, 43751 steps, wall=4481s, no NaN, CFL≤0.13, ibCG≤2 during
stroke, KE monotone-decaying in free decay): **every frame is clean** — no checkerboard at any
scale, including the blade-withdrawal window (frame 240) and the long free-decay tail.
R1+R2 alone (coarse 2h markers + Wendland-C2 volume kernel + reprojection) are the complete fix.
The resolved vortices are sharper and KE is higher than the AV-contaminated baseline.
Output: `ns_spoon_frames/`, `ns_spoon_flow_frames/`, `spoon.gif`, `spoon.png`.

## Notes / honest scope
- The lagged `ppe_div_damping` is kept (still load-bearing for stability); R1 reduces the
  divergence it has to mop up. R2 is the change that actually cleans the frames; R1 is
  correctness insurance.
- Marker-spacing/kernel translational-invariance *optimality* is a Cartesian result; on triangles
  these are smoothing/conditioning upgrades, not quantitative guarantees.
- M5 (fluid-side grid-scale dissipation) turned out to be a **non-issue** once M3+M4 were fixed:
  the checkerboard was never being *generated* by the flow itself; it was entirely an IB-transfer
  artifact. With the smooth kernel, the solver's natural (physical) viscosity at Re=500 is
  sufficient to keep the flow clean.
