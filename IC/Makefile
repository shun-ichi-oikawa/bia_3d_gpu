all: bbb_maker

FC=f95

bbb_maker: bbb_maker_REV.f90
	$(FC) -O3 -o bbb_maker bbb_maker_REV.f90

clean:
	rm bbb_maker
