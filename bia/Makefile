all: bia

FC = gfortran -O3
NV = nvcc -O3 -arch=sm_86

OBJS = find_dXYUV_kernel.o invariance.o sincos.o

bia: bia.f90 $(OBJS) Makefile
	$(FC) -o bia bia.f90 $(OBJS) -lcublas -lcudart -lstdc++

find_dXYUV_kernel.o: find_dXYUV_kernel.cu bia-rkf.cu Makefile
	$(NV) -c find_dXYUV_kernel.cu

invariance.o: invariance.cu Makefile
	$(NV) -c invariance.cu

sincos.o: sincos.c Makefile
	gcc -O3 -c sincos.c

clean:
	rm -f bia *.o *.mod
