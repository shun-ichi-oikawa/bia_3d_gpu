all: dim-cu

FC=f95 -O3 -lcublas
FC=ifort -O3 -xhost -parallel -lcublas
NV=nvcc -O3

GENCODE = \
-gencode=arch=compute_13,code=sm_13 -gencode=arch=compute_13,code=compute_13 \
-gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_20,code=compute_20

#dim-cu: dim_3D.f90 Invariance.f90 rkf_kernel.o Makefile
#	$(FC) -o dim-cu dim_3D.f90 rkf_kernel.o -lcudart -lcutil_x86_64 -lstdc++

dim-cu: dim_3D.f90 Invariance.f90 rkf_kernel.o fortran_thunking.o Makefile
	$(FC) -o dim-cu dim_3D.f90 rkf_kernel.o fortran_thunking.o -lcudart -lcutil_x86_64 -lstdc++

rkf_kernel.o: rkf_kernel.cu Makefile
	$(NV) $(GENCODE) -c rkf_kernel.cu

#fortran.o:/usr/local/cuda/src/fortran.c
#	$(NVC) /usr/local/cuda/src/fortran.c

fortran_thunking.o:	fortran_thunking.c
	$(NV) $(GENCODE) -c fortran_thunking.c -I/usr/local/cuda/src

clean:
	rm dim-cu rkf_kernel.o fortran_thunking.o
