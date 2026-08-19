# Lab4: AMSS-NCKU

Trimmed AMSS-NCKU numerical relativity lab for the single configuration in
`AMSS_NCKU_Input.py`:

- CPU evolution executable: `ABE`
- GPU evolution executable: `ABEGPU`
- Initial data: `Ansorg-TwoPuncture`
- Evolution equation: vacuum `BSSN`
- Grid: `Patch`, cell-centered, equatorial symmetry
- Finite difference: 4th order

Both CPU and GPU evolution paths are kept. Old verification scripts and
reference case files have been removed from this directory.

## Environment

The container is a Debian 13 box with the toolchain preinstalled. Two
architectures are shipped:

| Architecture | Hardware | Default toolchain | Executables |
| ------------ | -------- | ----------------- | ----------- |
| `linux/arm64` | Kunpeng 920B | GNU 14 + OpenMPI 5, `AMSS_ENABLE_GPU=OFF` | `TwoPunctureABE`, `ABE` |
| `linux/amd64` | x86_64 + NVIDIA A100 MIG | GNU 13 + OpenMPI 5 + CUDA 12.4, `sm_80` | `TwoPunctureABE`, `ABE`, `ABEGPU` |

x86/A100 uses **CUDA 12.4** (Debian 13 package). NVIDIA driver is **not**
in the image; the host injects it via NVIDIA Container Toolkit.

### Available compilers (advanced comparison)

The image ships several compiler/MPI combos beyond the default. Switch
via `CXX` / `FC` / `CUDACXX` / `MPI_CXX_COMPILER` env vars before
`compile.sh`.

- amd64: `icpx`/`ifx` (after `. /etc/profile.d/oneapi.sh`),
  `clang++`/`flang-19`, MPICH (`mpicxx.mpich` / `mpiexec.mpich`)
- arm64: `armclang++`/`armflang` (after `module load acfl/24.10.1`),
  MPICH (`mpicxx.mpich` / `mpiexec.mpich`)

CMake caches compilers, so switch toolchains with a fresh build dir.

### Switching MPI

Default is OpenMPI. To switch to MPICH (advanced comparison only):

```bash
export MPI_CXX_COMPILER=/usr/bin/mpicxx.mpich
export MPIEXEC_EXECUTABLE=/usr/bin/mpiexec.mpich
export AMSS_MPIEXEC=/usr/bin/mpiexec.mpich
rm -rf $AMSS_BUILD_DIR      # CMake cache is OpenMPI-specific; MUST use a fresh build dir
./compile.sh
```

Switching MPI **requires** a fresh `AMSS_BUILD_DIR` — the existing
`CMakeCache.txt` records the previous MPI's wrapper path and will fail
at configure time.

### CUDA-aware MPI

`AMSS_MPI_CUDA_AWARE` defaults to `0`. Debian OpenMPI/MPICH do not
advertise CUDA-aware support and the current AMSS code path uses
host-staging. Keep it at `0` for this checkout: the CUDA-aware interpolation
branches are incomplete and are not currently a supported build/test path.
If that code is repaired, pass `-DAMSS_MPI_CUDA_AWARE=1` to `compile.sh` only
after verifying that the MPI implementation actually supports device buffers.

### A100 MIG and `sm_80`

x86/A100 nodes expose MIG instances rather than an entire physical GPU. The
course partitions provide `1g.10gb` instances from A100 80GB GPUs and
`1g.5gb` instances from A100 40GB GPUs. MIG does not change the Ampere CUDA
architecture or compute capability, so the image's `nvcc` targets `sm_80`.
Use `nvidia-smi -L` to confirm the allocated MIG profile and `cuobjdump` to
inspect the architectures embedded in a GPU executable.

## Build

```bash
./compile.sh
```

This builds:

- `build/ABE`
- `build/ABEGPU` (only when `AMSS_ENABLE_GPU=ON`, which is the default on
  the amd64 image and `OFF` on the arm64 image)
- `build/TwoPunctureABE`

For a faster debug build:

```bash
./compile.sh -DAMSS_OPT='-O0'
```

## Run

```bash
./run.sh
```

`compile.sh` and `run.sh` resolve the repository from their own script path,
so an OJ runner may invoke them from a different working directory. The OJ
must run `compile.sh` before `run.sh`; `run.sh` intentionally does not compile
and checks that `ABE` and `TwoPunctureABE` already exist. With no extra CMake
flags, the default CPU build uses the validated POINTWISE/B=2 BSSN path and
TwoPuncture OpenMP support. The shipped input remains the 4.0-time smoke
case; a 40-step submission must change that setting in the OJ wrapper and
restore it afterward.

Optional TwoPuncture cache:

```bash
./run.sh --twop-cache
```

The run driver writes results under:

```text
GW250118/AMSS_NCKU_output/
GW250118/figure/
```

## Correctness check

```bash
./check.sh
```

`check.sh` resolves `RESULT_DIR` against `AMSS_OUTPUT_ROOT` (or the lab
root if unset), and `GOLDEN_DIR` against the lab root. The shipped
`golden/` directory is used by default. Pass an explicit `RESULT_DIR` to
check a non-default run directory. During a short DevPod/debug run the
checker may report `FINAL: PASS` on a matched prefix and print a warning when
the result has fewer timesteps than golden; that is not a complete formal
correctness result. Before recording a production PASS, confirm that the
matched timestep count equals the golden count. See
`python3 scripts/check_result.py --help` for details.

## Main Files

- `AMSS_NCKU_Input.py`: the fixed run parameters
- `AMSS_NCKU_Program.py`: run driver
- `scripts/setup.py`, `scripts/numerical_grid.py`,
  `scripts/generate_TwoPuncture_input.py`, `scripts/renew_puncture_parameter.py`:
  parfile generation
- `scripts/makefile_and_run.py`: launches `TwoPunctureABE` and the configured
  `ABE` or `ABEGPU` executable through MPI
- `scripts/plot_xiaoqu.py`, `scripts/plot_GW_strain_amplitude_xiaoqu.py`:
  post-run plots
- `scripts/check_result.py`: validates simulation output against golden results
- `src/`: source files still needed to compile `ABE`, `ABEGPU`, and `TwoPunctureABE`
