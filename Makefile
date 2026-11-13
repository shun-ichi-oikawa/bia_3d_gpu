all: bia

# --- Fortran Compiler Selection ---
# Standard GNU Fortran (Default for CUDA 11.8+ / modern Linux)
FC = gfortran -O3
# Intel Fortran (Used for legacy GTX 480 / CUDA 4.2 benchmark)
# FC = ifort -O3 -xhost -pc80

FS = bia.f90 Invariance_CPU.f90 Invariance_QPU.f90 \
find_dXYUV_RKF.f90 find_hyperbola.f90 find_parabola.f90 find_ellipse.f90

NV = nvcc -O3 -Xcompiler -march=native -Xcompiler -Wall

# --- GPU Architecture Selection ---
# Modern CUDA 11.x / 12.x (Adjust sm_XX to match your GPU, e.g., sm_75, sm_80, sm_86)
GENCODE = -arch=sm_75
# Legacy CUDA 4.2 (NVIDIA GeForce GTX 480 / Fermi)
# GENCODE = -arch=sm_20

# --- Linker Flags ---
# Modern CUDA (CUDA 5.0+)
LIB = -lcudart -lcublas -lstdc++
# Legacy CUDA 4.2 (Requires CUDA SDK cutil)
# LIB = -lcudart -lcutil_x86_64 -lcublas -lstdc++

bia: $(FS) find_dXYUV_kernel.o invariance.o sincos.o Makefile
	$(FC) -o bia bia.f90 find_dXYUV_kernel.o invariance.o sincos.o $(LIB)

sincos.o: sincos.c Makefile
	$(NV) $(GENCODE) -c sincos.c

invariance.o: invariance.cu Makefile
	$(NV) $(GENCODE) -c invariance.cu

find_dXYUV_kernel.o: find_dXYUV_kernel.cu bia-rkf.cu Makefile
	$(NV) $(GENCODE) -c find_dXYUV_kernel.cu

clean:
	rm -f bia *.o *.mod
