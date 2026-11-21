# bia_3d_gpu

A 3D hybrid Fortran/CUDA solver for gravitational and electrostatic $N$-body simulations using the Binary Interaction Approximation (`BIA`) alongside a direct $N$-body integration benchmark (`DIM`).

* **Version**: `v1.0.1` (Upgraded from legacy CUDA 4.2 to CUDA 11.8 with multi-architecture support `sm_61`–`sm_86`; legacy CUDA 4.2 release `v1.0.0` is archived on Zenodo at DOI `10.5281/zenodo.24131822`).

---

## Overview

`bia_3d_gpu` accelerates three-dimensional $N$-body simulations by decomposing pairwise interactions into analytical two-body Keplerian/Coulomb orbits (`ALG`: ellipse, parabola, hyperbola) or adaptive Runge–Kutta–Fehlberg (`RKF`) integrations on NVIDIA GPUs in strict `FP64` precision.

### Repository Structure
* **`bia/`**: Main Binary Interaction Approximation solver (`bia.f90`, `find_dXYUV_kernel.cu`, `bia-rkf.cu`, `invariance.cu`, `find_ellipse.f90`, `find_parabola.f90`, `find_hyperbola.f90`, `find_dXYUV4.f90`, `find_dXYUV_RKF.f90`, `Invariance_CPU.f90`, `Invariance_QPU.f90`, `sincos.c`, `Makefile`).
* **`dim/`**: Direct Integration Method (`DIM`) GPU reference solver (`dim_3D.f90`, `rkf_kernel.cu`, `Invariance.f90`, `fortran_thunking.c`, `Makefile`).
* **`IC/`**: Initial condition generator (`bbb_maker_REV.f90`, `Makefile`).
* **`OUT/`**: Output directory for initial conditions (`*.bbb`), direct integration baselines (`*.DIM`), and BIA results (`*.BIA`).

---

## Requirements & Building

* **GPU**: NVIDIA GPU with Compute Capability `sm_61` or higher (`Makefile` preconfigured for `sm_61`, `sm_70` [Volta / TITAN V], `sm_75`, and `sm_86`). Full-rate `FP64` hardware is recommended.
* **Software**: CUDA Toolkit 11.x (verified with CUDA 11.8), `gfortran` (`f95`), and `gcc`.

```bash
mkdir -p OUT
cd IC  && make && ./bbb_maker   # 1. Generate initial conditions in OUT/
cd ../dim && make && ./dim-cu   # 2. Run Direct Integration Method (DIM) baseline
cd ../bia && make && ./bia      # 3. Run Binary Interaction Approximation (BIA)
```

---

## Performance Benchmarks (N = 1,332 Particles, FP64)

Measured on **NVIDIA TITAN V** (Volta `sm_70`, 12 GB HBM2, CUDA 11.8):

| Solver / Mode | Kernel / Solver Time | Total Wall-Clock Time | Notes |
| :--- | :---: | :---: | :--- |
| **`dim-cu` (Direct Integration, GPU)** | — | `3.971 s` | Reference baseline ($\Delta E / E_0 = 2.44 \times 10^{-15}$) |
| **`bia` (`ALG` / `GPU` — Analytical)** | **`0.191 s`** | **`0.664 s`** | **~6.0x faster than DIM**, position diff $\mathcal{O}(10^{-13})$ |
| **`bia` (`ALG` / `CPU` — Analytical)** | `0.451 s` | `0.944 s` | Single-threaded Fortran CPU reference |
| **`bia` (`RKF` / `CPU` — Numerical)** | `1.951 s` | `2.423 s` | Adaptive Runge–Kutta–Fehlberg (CPU) |

*(Note: While the `RKF` / `GPU` mode was verified on legacy CUDA 4.2 / GTX 480 and CUDA 11.8 / RTX 3070, on Volta `sm_70` / TITAN V architectures the per-pair sequential kernel dispatch with independent thread scheduling incurs substantial latency; the analytical `ALG` / `GPU` mode is recommended for production runs.)*

---

## Author, Citation & License

* **Author**: Shun-ichi Oikawa (及川 俊一), Invited Teacher, Faculty of Engineering, Hokkaido University, Sapporo 060-8628, Japan
* **ORCID**: [https://orcid.org/0009-0001-5603-2787](https://orcid.org/0009-0001-5603-2787)
* **Contact**: `lopia312oikawa@gmail.com`
* **Zenodo Archive (v1.0.0)**: DOI `10.5281/zenodo.24131822`
* **License**: GNU General Public License v3.0 (GPL-3.0)
