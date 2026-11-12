# Binary Interaction Approximation (BIA) — 3D Multi-Body Hybrid Fortran/CUDA Code

[![DOI](https://zenodo.org)](https://doi.org)

A three-dimensional multi-body simulation package implementing the **Binary Interaction Approximation (BIA)** in a hybrid **Fortran 90 / NVIDIA CUDA C++** architecture.

## Key Features
* **Exact Momentum-Conserving Formulation:** Evaluates two-body relative displacements and velocity changes in the local orbital plane $(\xi, \eta, \zeta)$ and updates laboratory-frame states symmetrically, conserving total linear momentum to machine precision ($\sim 10^{-20}$).
* **Dual Solver Modes (`ALG` vs. `RKF`):**
  * `ALG` (*Exact Analytical Solver*): Integration-free exact Keplerian trajectory calculation (hyperbolic, parabolic, and elliptic) using Cardano's cubic formula and minimum-error-preserving Newton-Raphson iterations.
  * `RKF` (*6th-Order Runge-Kutta-Fehlberg*): Adaptive-step numerical ODE integration of the relative two-body equations of motion.
* **Hybrid CPU / GPU Execution (`PLATFM`):** Supports massively parallel execution on NVIDIA GPUs (`GPU`, CUDA C++) as well as multi-core host execution (`CPU`, Fortran).
* **Test-Particle Ensemble Mode (`N_end = 1`):** Efficiently computes $\Delta t = 1$ trajectories of a test particle placed at the origin across large initial-condition ensembles (`N_calc = 10**N_calc_EXP`) within the Debye sphere, directly yielding velocity-space diffusion and dynamical friction without small-angle Coulomb logarithm cutoffs.

## Physical Scaling & Specific Energy Convention
* **Dimensionless Interaction Parameter (`eps_th`):** Defined as $\varepsilon_{\text{th}} = \pm 2 \times 3 \times 10^{-7}$ (twice the single-particle coupling parameter due to reduced mass $\mu = m/2$, with negative sign for attractive interactions and positive sign for repulsive interactions).
* **Velocity-Squared Specific Energy (No $1/2$ Factor):** Throughout the entire code suite (`IC`, `dim`, and `bia`), kinetic and potential energies are normalized per unit mass with velocity-squared dimensions ($[v^2]$) **without the $1/2$ factor**:
  $$K = \sum_{i=1}^{N_p} |\mathbf{v}_i|^2, \qquad U = \sum_{i < j} \frac{2\varepsilon_{\text{th}}}{|\mathbf{r}_i - \mathbf{r}_j|}, \qquad E = K + U.$$

## Verified 4-Mode Benchmark ($N_p = 1,332$, $\Delta t = 1$)
Benchmarked on NVIDIA GeForce GTX 590 (CUDA 4.2) + Intel Xeon E5-2687W (`ifort -O3 -xhost -pc80`) against full GPU Direct Integration (`dim-cu`, total runtime `8.402 s`, relative energy error $-9.99 \times 10^{-16}$):

| Mode (`PLATFM` $\times$ `METHOD`) | Pairwise Compute Time | Total Wall Time | Position RMS Error (`std_r`) | Velocity RMS Error (`std_v`) |
| :--- | :---: | :---: | :---: | :---: |
| **`GPU` $\times$ `ALG`** | **0.342 s** | **1.028 s** | $1.28 \times 10^{-12}$ | $3.08 \times 10^{-12}$ |
| **`CPU` $\times$ `ALG`** | 0.528 s | 1.170 s | $1.56 \times 10^{-12}$ | $3.44 \times 10^{-12}$ |
| **`GPU` $\times$ `RKF`** | 2.255 s | 2.899 s | $7.47 \times 10^{-13}$ | $1.32 \times 10^{-12}$ |
| **`CPU` $\times$ `RKF`** | 3.012 s | 3.661 s | $7.47 \times 10^{-13}$ | $1.32 \times 10^{-12}$ |

## Repository Structure & End-to-End Workflow

### Directory Layout
* `IC/` : Initial condition generator (`bbb_maker_REV.f90`, `Makefile`).
* `dim/` : Direct Integration Method reference solver (`dim_3D.f90`, `rkf_kernel.cu`, `fortran_thunking.c`, `Invariance.f90`, `Makefile`).
* `OUT/` : Directory for binary initial conditions (`.bbb`), exact references (`.DIM`), and BIA outputs (`.BIA`, `BIA_10^k.*`).
* `./` (Root) : Main 3D BIA hybrid Fortran/CUDA solver (`bia.f90`, `find_dXYUV_kernel.cu`, `bia-rkf.cu`, analytical Kepler solvers, invariance auditors, and `Makefile`).

### Step 1: Generate Initial Conditions (`IC/`)
Configure `N_L` (shell layers, $N_p = (2N_L + 1)^3 + 1$) and `N_calc_EXP` ($10^{\text{N\_calc\_EXP}}$ ensembles) in `IC/bbb_maker_REV.f90`, then build and run:
```bash
cd IC
make && ./bbb_maker
```
This generates `../OUT/0001332_10^0.bbb` (for `N_L = 5, N_calc_EXP = 0`) with the test particle at the origin and field particles ordered outward by shell.

### Step 2: Generate High-Precision Reference (`dim/`, Optional)
For full-system verification (`N_end == N_p`), configure `NP` in `dim/dim_3D.f90` and run the GPU 6th-order RKF solver:
```bash
cd ../dim
make && ./dim-cu
```
This writes the ground-truth reference trajectory to `../OUT/0001332_10^0.DIM`.
## Repository Structure & Quick Start (`IC` $\to$ `dim` $\to$ `bia`)

* **`./` (Root):** 3D Multi-Body Hybrid Fortran/CUDA BIA Solver (`bia.f90`, `find_dXYUV_kernel.cu`, `Makefile`, etc.)
* **`IC/`:** Layered-shell initial-condition binary generator (`bbb_maker_REV.f90`, `Makefile`)
* **`dim/`:** GPU-accelerated Direct Integration Method (RKF) reference solver (`dim_3D.f90`, `rkf_kernel.cu`, `Makefile`)
* **`OUT/`:** Output directory for `.bbb` (initial), `.DIM` (reference), and `.BIA` / `BIA_*` (results)

### Step 1: Generate Initial Conditions (`IC/`)
Set `N_L = 5` ($N_p = 1,332$) and `N_calc_EXP = 0` in `IC/bbb_maker_REV.f90`:
```bash
cd IC && make && ./bbb_maker && cd ..
```

### Step 2: Generate Reference Solution (`dim/`, Optional)
Required only when running full-system accuracy comparisons (`N_end = N_p`):
```bash
cd dim && make && ./dim-cu && cd ..
```

### Step 3: Run 3D Multi-Body BIA Solver (`./bia`)
Set `METHOD` (`'ALG'` or `'RKF'`), `PLATFM` (`'GPU'` or `'CPU'`), `N_p = 1332`, and `N_end = N_p` (or `N_end = 1` for test-particle ensembles) in `bia.f90`:
```bash
make clean && make && ./bia

## Citation & License
* **This Software Archive (3D Hybrid Fortran/CUDA, Zenodo DOI):** Shun-ichi Oikawa, *Binary Interaction Approximation (BIA) — 3D Multi-Body Hybrid Fortran/CUDA Code (v1.0.0)*, Zenodo (2026). DOI: [10.5281/zenodo.24131822](https://doi.org/10.5281/zenodo.24131822)
* **Sequential C/GSL Reference Archive (Zenodo DOI):** Shun-ichi Oikawa, *Binary Interaction Approximation (BIA) Code Package (v1.0.1)*, Zenodo (2026). DOI: [10.5281/zenodo.24122562](https://doi.org)
* **Foundational Papers:** S. Oikawa et al., *Plasma and Fusion Research* (2011, 2013, 2014).
* **License:** GNU General Public License v3.0 (GPL-3.0).
