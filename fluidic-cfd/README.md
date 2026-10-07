# 4×4 fluidic matrix multiplier — full 3D CFD

Extract the entire ZIP before running any scripts. To confirm that every
packaged file is present and unmodified, run `python3 verify-package.py` from
the extracted `fluidic-cfd` directory. It checks `MANIFEST.sha256` without any
third-party dependencies. A file edited or regenerated after extraction will
naturally differ from its original packaged checksum.

For the interactive results, open `fluidic-matrix-multiplier-3d.html` in a
modern browser. To rebuild or run a new experiment, follow **Reproduce or modify**
below. The archive contains the complete source and data; a platform-specific
compiled executable is built locally by `run.py`.

This project solves fluid flow throughout a three-dimensional channel network.
It does **not** use a hydraulic-resistor network, prescribed parabolic channel
profiles, or animated particles as a substitute for a flow solution.

Four upper pressure manifolds connect through sixteen square vertical channels
to four lower output manifolds. A run produces a three-component velocity field
and a pressure field in every fluid voxel, including all junctions. The first
device has **nonnegative geometric weights**. Inputs may have either sign as
gauge pressures. Negative matrix coefficients would require an additional
differential-output architecture; that is not implemented here.

## What was solved

The intended continuum model is the isothermal, Newtonian Navier–Stokes system:

\[
\partial_t\rho+\nabla\cdot(\rho\mathbf u)=0,\qquad
\partial_t(\rho\mathbf u)+\nabla\cdot(\rho\mathbf u\mathbf u)
=-\nabla p+\nabla\cdot\boldsymbol\tau.
\]

A D3Q19 BGK lattice-Boltzmann discretization evolves the distributions in the
full volume, including the quadratic equilibrium terms that retain convective
inertia. The pressure relation is \(p=c_s^2\rho\), with \(c_s^2=1/3\) in lattice
units. This is a **weakly compressible, low-Mach approximation** to incompressible
water flow, not an exact incompressible projection or a compressible acoustic
model of real water. Low lattice Mach numbers are checked, not assumed.

All solid surfaces use halfway bounce-back no-slip walls. The fluid has
dynamic viscosity 0.001 Pa·s and reference density 1,000 kg/m³. Walls are rigid.
The model excludes thermal transport, bubbles, cavitation, surface chemistry,
elastic walls, pump dynamics, measurement noise, and fabrication tolerances.

The pressure ports use nonequilibrium extrapolation from the adjacent interior
plane. Input reservoir pressures are prescribed; all four output pressures
are 0 Pa gauge. Boundary velocity is extrapolated rather than prescribed.
There is no lumped compliance, imposed channel flow rate, or active feedback.
The ideal pressure reservoirs may either supply or absorb fluid; in particular,
a grounded input can receive backflow through the connected network.

## Geometry and operating point

- Fluid domain envelope: 880 × 880 × 480 µm.
- Four upper input manifolds: 120 × 120 µm cross-section, 800 µm length.
- Four lower output manifolds: the same section and length, perpendicular to the inputs.
- Upper/lower manifold center heights: 360 / 100 µm.
- Sixteen vertical square channels overlap the manifolds; their center-to-center
  length is 260 µm, and the unobstructed gap between manifolds is 140 µm.
- Connecting widths are specified in `widths.txt`, in units of 10 µm.
- Mixed input pressures: [0.8, 0.4, 0.2, 0.6] Pa.
- Signed test pressures: [0.8, 0.4, −0.2, 0.6] Pa.

The channels and sharp junctions are grid-aligned rectangular solids. Geometry
is identical at 10 and 5 µm spacing: the finer grid subdivides each coarse cell
into eight. The solver includes 100,112 and 800,896 fluid voxels respectively.
The four smallest-channel cells across the coarse grid are deliberately checked
against eight across the fine grid. This is a prototype CFD calculation, not a
fabrication-ready or certified mesh-independent design.

## The matrix being multiplied

Four independent 1 Pa basis-pressure runs calibrate the device's actual transfer
matrix \(G\). Column \(i\) is the vector of four output flows when only input
\(i\) has nonzero pressure. Its units are nL/s per Pa. Separate combined-input
runs test \(\mathbf Q\approx G\mathbf P\).

The full 3D solution determines these coefficients. The channel widths are
geometric parameters, **not claimed to be prescribed numerical weights**.
Manifold pressure losses, junction flow, and interactions between branches all
contribute to the calibrated matrix. This prototype characterizes the matrix
of one device; it does not inverse-design arbitrary requested coefficients.

The displayed flows are **mass-equivalent volumetric rates**: mass flux divided
by the reference density of water. For conservation, the solver sums every
lattice link crossing an interior port plane, including diagonal links.
This avoids an inconsistent center-velocity quadrature. Small differences
from local volumetric flux are possible in a weakly compressible method.

## Verification

Exact measured results are in `results/validation.json`. The checks include:

1. All device runs reach a relative L2 velocity-change tolerance of 10⁻⁶ for
   three consecutive checks, spaced 0.002 s apart. The zero-input case uses an
   absolute lattice speed floor of 10⁻¹² to avoid dividing by roundoff.
2. Flux conservation across all eight ports is checked independently of the
   velocity-change criterion.
3. Four basis runs are compared with independently solved mixed and signed
   pressure runs. No linear combination of fields is substituted for those runs.
4. The mixed-input case is solved on both 10 and 5 µm grids. The reported change
   is **grid sensitivity**, not a continuum error bound; two grids do not prove
   an asymptotic convergence order.
5. An independent straight square duct is checked against the analytical
   fully developed flow series. Its pressure gradient is measured from the
   interior CFD pressure field, and its volumetric flow is integrated on the
   middle cross-section. Both resolutions are tested.
6. Zero input is checked for spurious flow. Maximum lattice Mach number is recorded.

The standalone viewer selects actual computed cases. Changing its case selector
does not run a new CFD solve. Its point cloud is a regular sample of volume cells;
the full data are included separately. Arrows represent computed local velocity
directions and a scaled magnitude, not particle trajectories. The transparent
channel outlines show the design dimensions, not a fitted surface mesh.

## Reproduce or modify

Requirements: Python 3 with NumPy, and `g++` with OpenMP. No commercial solver is
required. Optional VTK Python bindings export a ParaView volume file. The solver
uses double precision; exported volume samples use 32-bit floats.

```sh
python3 run.py --case all --threads 8
```

This rebuilds and runs all basis, combined-input, signed-input, zero-input,
resolution, and duct checks. It then rebuilds the viewer. Runtime depends on CPU
and memory bandwidth; a fine-grid run is substantially more expensive than a
coarse run. Exit status 5 means that the maximum iteration count was reached
without convergence; the result must not be treated as converged.

For a new experiment, edit the four pressures in `inputs.txt` (normalized to
the pressure scale, each between −1 and 1), and optionally the sixteen widths
in `widths.txt`. Widths must be even integers from 0 through 10; 0 blocks a
connecting channel. Each unit is 10 µm.

```sh
python3 run.py --case custom --scale 2 --pressure 1 --threads 8
python3 export-vtk.py results/custom-s2
```

Changing geometry requires recalibrating all four basis runs before reusing the
matrix. The current viewer is built for the bundled cases; a custom experiment's
full field is inspected in ParaView or by reading its binary arrays. For a true
4×4 times 4×4 product, run four input pressure vectors as the columns of the
second matrix; the four resulting output vectors are the product columns,
within the measured linearity and discretization errors.

## Files

- `solver.cpp`: full volumetric LBM solver, boundary conditions, convergence
  diagnostics, conservative port-flux integration, and binary field export.
- `run.py`: reproducible build/run driver.
- `analyze.py`: independent matrix, grid, duct, and conservation checks.
- `results/*-field.bin`: one seven-float little-endian record per fluid voxel:
  x µm, y µm, z µm, pressure Pa, ux mm/s, uy mm/s, uz mm/s.
- `results/*.json`: solver metadata, settings, diagnostics, and port flows.
- `results/*-history.csv`: iteration histories. These are solver convergence
  histories, not predictions of real-water acoustic transients.
- `results/mix-s2.vti`: complete fine-grid field for ParaView. Threshold
  `fluid_mask` at 1, then color by pressure or velocity. Coordinates are µm.
- `viewer-template.html`, `viewer.js`, `build-viewer.py`: interactive viewer sources.
- `vendor/three.min.js`: Three.js 0.160.0, MIT licensed.

Both viewers bundle Three.js 0.160.0 (retrieved from its version-pinned CDN) and
the computed data. No network access is required. A modern browser with WebGL
and the standard DecompressionStream API is required.

## Method references

- Zou & He, *On pressure and velocity boundary conditions for the lattice
  Boltzmann BGK model*, Physics of Fluids 9, 1591 (1997).
  https://doi.org/10.1063/1.869307
- Guo, Zheng & Shi, *Non-equilibrium extrapolation method for velocity and
  pressure boundary conditions in the lattice Boltzmann method*, Chinese Physics
  11, 366 (2002). https://doi.org/10.1088/1009-1963/11/4/310
- Nash et al., *Choice of boundary condition for lattice-Boltzmann simulation
  of moderate-Reynolds-number flow in complex domains*, Physical Review E 89,
  023303 (2014). https://doi.org/10.1103/PhysRevE.89.023303

The implementation and tests here are new project code. These references
describe the method family; they do not independently validate this device.
