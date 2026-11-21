all: dim-cu

FC=f95 -O3 -lcublas
#C=ifort -O3 -xhost -parallel -lcublas
NV=nvcc -O3

GENCODE = \
-gencode=arch=compute_86,code=sm_86 -gencode=arch=compute_86,code=compute_86

#dim-cu: dim_3D.f90 Invariance.f90 rkf_kernel.o Makefile
#	$(FC) -o dim-cu dim_3D.f90 rkf_kernel.o -lcudart -lstdc++

dim-cu: dim_3D.f90 Invariance.f90 rkf_kernel.o fortran_thunking.o Makefile
	$(FC) -o dim-cu dim_3D.f90 rkf_kernel.o fortran_thunking.o -lcudart -lstdc++

rkf_kernel.o: rkf_kernel.cu Makefile
	$(NV) $(GENCODE) -c rkf_kernel.cu

#fortran.o:/usr/local/cuda/src/fortran.c
#	$(NVC) /usr/local/cuda/src/fortran.c

fortran_thunking.o:  /usr/local/cuda/src/fortran_thunking.c
	$(NV) $(GENCODE) -DCUBLAS_GFORTRAN -c /usr/local/cuda/src/fortran_thunking.c -I/usr/local/cuda/src

clean:
	rm dim-cu rkf_kernel.o fortran_thunking.o
