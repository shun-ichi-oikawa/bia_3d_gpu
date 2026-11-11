#all: bia

#C = f95 -O3
#C = ifort -O3 -xhost -parallel
FC = ifort -O3 -xhost -pc80

FS = bia.f90 Invariance_CPU.f90 Invariance_QPU.f90 \
find_dXYUV_RKF.f90 find_hyperbola.f90 find_parabola.f90 find_ellipse.f90

NV = nvcc -O3 -Xcompiler -march=native -Xcompiler -Wall

LIB = -lcudart -lcutil_x86_64 -lcublas -lstdc++

GENCODE = -arch=sm_20

#-gencode=arch=compute_13,code=sm_13 -gencode=arch=compute_13,code=compute_13 \
#-gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_20,code=compute_20

bia: $(FS) find_dXYUV_kernel.o invariance.o sincos.o Makefile
	$(FC) -o bia bia.f90 find_dXYUV_kernel.o invariance.o sincos.o $(LIB)

sincos.o: sincos.c Makefile
	$(NV) $(GENCODE) -c sincos.c

invariance.o: invariance.cu Makefile
	$(NV) $(GENCODE) -c invariance.cu

find_dXYUV_kernel.o: find_dXYUV_kernel.cu bia-rkf.cu Makefile
	$(NV) $(GENCODE) -c find_dXYUV_kernel.cu

clean: 
	rm bia *.o
