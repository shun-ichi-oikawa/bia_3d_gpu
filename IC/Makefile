all: bbb_maker

# Standard GNU Fortran (Default for portability)
FC = gfortran

bbb_maker: bbb_maker_REV.f90
	$(FC) -O3 -o bbb_maker bbb_maker_REV.f90

clean:
	rm -f bbb_maker *.mod
