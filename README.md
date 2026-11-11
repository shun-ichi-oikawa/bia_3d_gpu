# Binary Interaction Approximation (BIA) — 3D Multi-Body Hybrid Fortran/CUDA Code

[![DOI](https://zenodo.org)](https://doi.org/10.5281/zenodo.24131822)

## Overview & Architecture
This repository contains the official 3D multi-body parallel implementation of the Binary Interaction Approximation (BIA) framework for high-precision N-body Coulombic plasma and gravitational simulations. Built as a **strictly momentum-conserving** hybrid Fortran 90 and CUDA C++ solver, it decomposes 3D multi-body dynamics into exact, pairwise 2D Keplerian orbital trajectories (`hyperbola`, `parabola`, and `ellipse`) without time-step discretization errors.

### Key Features
* **Strictly Momentum-Conserving Formulation:** Enforces exact pairwise action-reaction symmetry across all interacting pairs, preserving total linear momentum down to $\mathcal{O}(10^{-21})$–$\mathcal{O}(10^{-19})$ without artificial center-of-mass drift.
* **Dual Solver Modes (`ALG` vs. `RKF`):** Features an integration-free exact analytical mode (`ALG`) solving Keplerian trajectories via Newton-Raphson iteration, alongside a numerical 6th-order Runge-Kutta-Fehlberg mode (`RKF`) for benchmarking.
* **Hybrid CPU/GPU Execution (`PLATFM = 'GPU'` / `'CPU'`):** Offloads $\mathcal{O}(N^2)$ pairwise relative orbital evaluations onto NVIDIA CUDA threads while utilizing Intel Fortran for host orchestration and precision diagnostics.
* **Triple-Precision Conservation Auditing:** Cross-verifies total linear momentum, angular momentum, and total energy ($KE + PE$) across GPU double precision (`CUDA`), CPU double precision (`DBLE`), and CPU quadruple precision (`QUAD`).

---

## Physical Scaling & `eps_th` Parameter
* **`eps_th` (Dimensionless Interaction Parameter):** Represents the Coulombic (or gravitational) interaction strength normalized by the particle mass, setting the characteristic relative thermal speed to $1$ and the mean inter-particle distance to $1$. Under typical fusion-grade proton plasma conditions, this normalized single-particle interaction parameter is approximately $3 \times 10^{-7}$.
* **Factor of $2$ (Reduced Mass Effect):** Solving the two-body relative equation of motion ($\mathbf{r}_{ij} = \mathbf{r}_i - \mathbf{r}_j$) divides by the reduced mass $\mu = m/2$, introducing an exact factor of $2$ into the relative acceleration equation ($\varepsilon_{\text{th}} = \pm 2 \times 3 \times 10^{-7} = \pm 6.0 \times 10^{-7}$).
* **Sign Convention:** Negative (`-6.0e-7`, i.e., `-2 * 3d-7`) denotes **attractive** interactions (gravitational systems or opposite-charge Coulombic attraction), whereas positive (`+6.0e-7`, i.e., `+2 * 3d-7`) denotes **repulsive** interactions (same-species fusion protons). This normalization enables exact analytical integration across macroscopic time steps of $\Delta t = 1$.

---

## Verified 4-Mode Benchmark Performance ($N_p = 1,332$)
Validated on a CUDA Fermi-architecture workstation against a Direct Integration Method (`DIM`) reference over a macroscopic step $\Delta t = 1$:

| Execution Mode (`PLATFM` $\times$ `METHOD`) | Pairwise Time | Total Wall Time | Position RMS Error (`std_r`) | Velocity RMS Error (`std_v`) | Momentum Drift (`sum(dU)`) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`GPU` $\times$ `ALG` (Integration-Free Exact)** | **`0.342 sec`** | **`1.028 sec`** | `1.28e-12` | `3.08e-12` | `7.71e-20` |
| **`CPU` $\times$ `ALG` (Integration-Free Exact)** | `0.528 sec` | `1.170 sec` | `1.56e-12` | `3.08e-12` | `-5.08e-21` |
| **`GPU` $\times$ `RKF` (6th-Order Numerical)** | `2.255 sec` | `2.899 sec` | `7.47e-13` | `3.08e-12` | `2.38e-19` |
| **`CPU` $\times$ `RKF` (6th-Order Numerical)** | `3.012 sec` | `3.661 sec` | `7.47e-13` | `3.08e-12` | `7.37e-20` |

* **Speedup:** The integration-free CUDA solver (`GPU` $\times$ `ALG`) achieves an **8.8$\times$ speedup** in pairwise evaluation over the traditional numerical integrator (`CPU` $\times$ `RKF`), while matching the `DIM` exact reference to **12 significant digits**.

---

## Repository Structure & Usage
* Core files include `bia.f90`, CUDA device kernels (`find_dXYUV_kernel.cu`, `bia-rkf.cu`), Fortran Keplerian solvers (`find_hyperbola.f90`, `find_parabola.f90`, `find_ellipse.f90`, `find_dXYUV4.f90`, `find_dXYUV_RKF.f90`), C math wrapper (`sincos.c`), and quadruple-precision auditors (`Invariance_CPU.f90`, `Invariance_QPU.f90`, `invariance.cu`).
* **Build & Run:** Configure `METHOD`, `PLATFM`, `N_p`, and `eps_th` in `bia.f90`, compile with `make clean && make`, and execute `./bia` using initial conditions in `../OUT/`.

---

## Citation & License
* **This Software Archive (3D Hybrid Fortran/CUDA, Zenodo DOI):** Shun-ichi Oikawa, *Binary Interaction Approximation (BIA) — 3D Multi-Body Hybrid Fortran/CUDA Code (v1.0.0)*, Zenodo (2026). DOI: [10.5281/zenodo.24131822](https://doi.org/10.5281/zenodo.24131822)
* **Sequential C/GSL Reference Archive (Zenodo DOI):** Shun-ichi Oikawa, *Binary Interaction Approximation (BIA) Code Package (v1.0.1)*, Zenodo (2026). DOI: [10.5281/zenodo.24122562](https://doi.org)
* **Foundational Papers:** S. Oikawa et al., *Plasma and Fusion Research* (2011, 2013, 2014).
* **License:** GNU General Public License v3.0 (GPL-3.0).
