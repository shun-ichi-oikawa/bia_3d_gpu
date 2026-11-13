all: dim-cu

# --- Fortran Compiler Selection ---
# Standard GNU Fortran (Default for portability)
FC = gfortran -O3 -lcublas
# Intel Fortran (Uncomment if using ifort on legacy/HPC systems)
# FC = ifort -O3 -xhost -parallel -lcublas

NV = nvcc -O3

# --- GPU Architecture Selection ---
# Modern CUDA 11.x / 12.x (Adjust sm_XX to match your GPU, e.g., sm_75, sm_80, sm_86)
GENCODE = -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_75,code=compute_75

# Legacy CUDA 4.2 (NVIDIA GeForce GTX 480 / Fermi)
# GENCODE = \
# -gencode=arch=compute_13,code=sm_13 -gencode=arch=compute_13,code=compute_13 \
# -gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_20,code=compute_20

# --- Linker Flags ---
# Modern CUDA (CUDA 5.0+)
LDFLAGS = -lcudart -lstdc++
# Legacy CUDA 4.2 (Requires CUDA SDK cutil)
# LDFLAGS = -lcudart -lcutil_x86_64 -lstdc++

dim-cu: dim_3D.f90 Invariance.f90 rkf_kernel.o fortran_thunking.o Makefile
	$(FC) -o dim-cu dim_3D.f90 Invariance.f90 rkf_kernel.o fortran_thunking.o $(LDFLAGS)

rkf_kernel.o: rkf_kernel.cu Makefile
	$(NV) $(GENCODE) -c rkf_kernel.cu

fortran_thunking.o: fortran_thunking.c
	$(NV) $(GENCODE) -c fortran_thunking.c -I/usr/local/cuda/src

clean:
	rm -f dim-cu rkf_kernel.o fortran_thunking.o *.mod
