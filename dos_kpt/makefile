all: dos_kpt.f90
	gfortran dos_kpt.f90 -o exe.x -ffree-line-length-none
#-stand f95
dft: dos_kpt_fixed.f90
	ifort dos_kpt_fixed.f90 -o fixed.x -O3

clean:
	rm *.x
