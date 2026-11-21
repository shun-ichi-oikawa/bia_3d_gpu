# bia_3d_gpu

A 3D hybrid Fortran/CUDA solver for gravitational and electrostatic $N$-body simulations using the Binary Interaction Approximation (`BIA`) alongside a direct $N$-body integration benchmark (`DIM`).

* **Version**: `v1.0.1` (Upgraded from legacy CUDA 4.2 to CUDA 11.8 with multi-architecture support `sm_61`–`sm_86`; legacy CUDA 4.2 release `v1.0.0` is archived on Zenodo at DOI `10.5281/zenodo.24131822`).

---

## Overview

`bia_3d_gpu` accelerates three-dimensional $N$-body simulations by decomposing pairwise interactions into analytical two-body Keplerian/Coulomb orbits (`ALG`: ellipse, parabola, hyperbola) or adaptive Runge–Kutta–Fehlberg (`RKF`) integrations on NVIDIA GPUs in strict `FP64` precision.

### Repository Structure
* **`bia/`**: Main Binary Interaction Approximation solver (`bia.f90`, `find_dXYUV_kernel.cu`, `bia-rkf.cu`, `invariance.cu`, `find_ellipse.f90`, `find_parabola.f90`, `find_hyperbola.f90`, `find_dXYUV4.f90`, `find_dXYUV_RKF.f90`, `Invariance_CPU.f90`, `Invariance_QPU.f90`, `sincos.c`, `Makefile`).
* **`dim/`**: Direct Integration Method (`DIM`) GPU reference solver using a 6-stage, 5th-order Runge–Kutta–Fehlberg (`RKF`) method (`dim_3D.f90`, `rkf_kernel.cu`, `Invariance.f90`, `fortran_thunking.c`, `Makefile`).
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
| **`dim-cu` (Direct Integration, GPU — 6-stage, 5th-order RKF)** | — | `2.946 s` | Reference baseline ($\Delta E / E_0 = 1.33 \times 10^{-15}$) |
| **`bia` (`ALG` / `GPU` — Analytical)** | **`0.193 s`** | **`0.680 s`** | **~4.3x faster than DIM**, position diff $\mathcal{O}(10^{-13})$ |
| **`bia` (`ALG` / `CPU` — Analytical)** | `0.445 s` | `0.917 s` | Single-threaded Fortran CPU reference, position diff $\mathcal{O}(10^{-13})$ |
| **`bia` (`RKF` / `GPU` — Numerical)** | `0.293 s` | `0.766 s` | Adaptive Runge–Kutta–Fehlberg (GPU) |
| **`bia` (`RKF` / `CPU` — Numerical)** | `1.934 s` | `2.413 s` | Adaptive Runge–Kutta–Fehlberg (CPU), position diff $\mathcal{O}(10^{-13})$ |

*(Note: Total wall-clock time for `bia` includes initial condition loading and multi-precision invariance verification across GPU, CPU, and Quadruple-precision routines; for production runs, the analytical `ALG` / `GPU` mode offers the highest throughput and machine-epsilon agreement with `DIM`.)*

---

## Author, Citation & License

* **Author**: Shun-ichi Oikawa (及川 俊一), Invited Teacher, Faculty of Engineering, Hokkaido University, Sapporo 060-8628, Japan
* **ORCID**: [https://orcid.org](https://orcid.org)
* **Contact**: `lopia312oikawa@gmail.com`
* **Zenodo Archive (v1.0.1 — Latest)**: DOI `10.5281/zenodo.24414593`
* **Zenodo Archive (v1.0.0 — Legacy CUDA 4.2)**: DOI `10.5281/zenodo.24131822`
* **License**: GNU General Public License v3.0 (GPL-3.0)


