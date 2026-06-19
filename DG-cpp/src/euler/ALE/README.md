# Body-fitted ALE DG Examples

This folder contains body-fitted ALE prototypes for compressible-Euler DG on
moving fitted meshes.  The code lives under `src/euler/ALE` because it reuses the
Euler DG solver and Euler AMR scene infrastructure.

The code is split by responsibility:

- `core/`: shared ALE infrastructure, including DG residual/update helpers,
  DG state interpolation, checkpoint/output utilities, body-fitted mesh
  generation, elastic solid elements, and reusable piston model code.
- `cases/Convergence.h/.cpp`: moving-mesh isentropic-vortex convergence test.
- `cases/Piston.h/.cpp`: prescribed moving-wall piston case.
- `cases/FSI.h/.cpp`: two-way mass-spring piston FSI case.
- `cases/Panel.h/.cpp`: two-way flexible top-wall modal panel case.
- `cases/TurekFlag.h/.cpp`: cylinder with an attached flexible flag benchmark
  surrogate on a body-fitted SDF mesh.
- `cases/BlastRod/`: high-pressure gas flow over a bottom-clamped upright
  elastic rod in a rectangular channel.

Each case `.cpp` contains its own `main` function; there is no shared dispatcher
and no separate `*_main.cpp` entry file.

All generated artifacts, including frames, still images, videos, and CSV
diagnostics, must be written under `out/`.

Build:

```sh
cmake -S DG-cpp -B DG-cpp/build
cmake --build DG-cpp/build --target \
  ale_convergence piston fsi panel \
  turek_flag blast_rod -j
```

Run the moving-mesh convergence check:

```sh
./DG-cpp/build/ale_convergence
```

Run a quick prescribed piston movie:

```sh
cd DG-cpp
./build/piston --quick
ffmpeg -y -framerate 25 -i out/piston_quick_frames/frame_%05d.ppm \
  -c:v libx264 -pix_fmt yuv420p -crf 16 out/piston_quick.mp4
```

Run the two-way mass-spring piston FSI case:

```sh
cd DG-cpp
./build/fsi --quick
ffmpeg -y -framerate 25 -i out/fsi_quick_frames/frame_%05d.ppm \
  -c:v libx264 -pix_fmt yuv420p -crf 16 out/fsi_quick.mp4
```

Run the two-way flexible-panel case:

```sh
cd DG-cpp
./build/panel --quick
ffmpeg -y -framerate 25 -i out/panel_quick_frames/frame_%05d.ppm \
  -c:v libx264 -pix_fmt yuv420p -crf 16 out/panel_quick.mp4
```

Run the cylinder-attached-flag case:

```sh
cd DG-cpp
./build/turek_flag --quick
ffmpeg -y -framerate 25 -i out/turek_flag_quick_frames/frame_%05d.ppm \
  -c:v libx264 -pix_fmt yuv420p -crf 16 out/turek_flag_quick.mp4
ffmpeg -y -i out/turek_flag_quick.ppm -frames:v 1 -update 1 out/turek_flag_quick.png
```

The FSI, panel, Turek-style flag, and blast-rod runs write CSV diagnostics under `out/`.

Run the high-pressure gas / upright elastic rod case. This case uses ALE motion
on a static, body-fitted SDF mesh with dense local resolution around the rod;
it does not perform dynamic AMR during the time loop.

```sh
cd DG-cpp
./build/blast_rod --quick
ffmpeg -y -framerate 25 -i out/blast_rod_quick_frames/frame_%05d.ppm \
  -c:v libx264 -pix_fmt yuv420p -crf 16 out/blast_rod_quick.mp4
ffmpeg -y -framerate 25 -i out/blast_rod_quick_nomesh_frames/frame_%05d.ppm \
  -c:v libx264 -pix_fmt yuv420p -crf 16 out/blast_rod_quick_nomesh.mp4
ffmpeg -y -i out/blast_rod_quick.ppm -frames:v 1 -update 1 out/blast_rod_quick.png
ffmpeg -y -i out/blast_rod_quick_nomesh.ppm -frames:v 1 -update 1 out/blast_rod_quick_nomesh.png
```

Run the dense 60-fps presentation version:

```sh
cd DG-cpp
./build/blast_rod
ffmpeg -y -framerate 60 -i out/blast_rod_frames/frame_%05d.ppm \
  -c:v libx264 -pix_fmt yuv420p -crf 16 out/blast_rod.mp4
ffmpeg -y -framerate 60 -i out/blast_rod_nomesh_frames/frame_%05d.ppm \
  -c:v libx264 -pix_fmt yuv420p -crf 16 out/blast_rod_nomesh.mp4
ffmpeg -y -i out/blast_rod.ppm -frames:v 1 -update 1 out/blast_rod.png
ffmpeg -y -i out/blast_rod_nomesh.ppm -frames:v 1 -update 1 out/blast_rod_nomesh.png
```
