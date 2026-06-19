# Common Infrastructure

`MeshGen.h` / `MeshGen.cpp` provide the shared unstructured triangular mesh
generator used by Navier-Stokes, Euler, and ALE-AMR drivers.

The generic entry point is `generateDistanceMesh`:

- `DistanceMeshSpec::signedDistance` defines the fluid domain; values are
  negative inside, zero on the boundary, and positive outside.
- `DistanceMeshSpec::targetSize` defines the local desired edge length.
- `DistanceMeshSpec::seedH` optionally sets a finer initial candidate-point
  spacing when `targetSize` has strongly localized small values.
- `DistanceMeshSpec::fixedPoints` pins exact boundary samples such as corners,
  circle rings, flag roots, or other geometry anchors.
- `DistanceMeshSpec::projectInside` is optional; without it, escaped nodes are
  projected back by a numerical signed-distance gradient.

The older Navier-Stokes helpers remain available as compatibility wrappers:

- `generateCylinderMesh` / `classifyEdges`
- `generateBowlMesh` / `classifyBowlEdges`

New body-fitted ALE-AMR benchmark geometries should build a `DistanceMeshSpec`
instead of adding solver-local mesh generators.
