dos: dos_traj.f90
	ifort dos_traj.f90 -o exe.x -O3 -static -stand f95
pos: xyz_to_POSCAR.f90
	ifort xyz_to_POSCAR.f90 -o xyz_to_POSCAR.x -O3 -static -stand f95
geo: xyz_to_GEOMETRY.f90
	ifort xyz_to_GEOMETRY.f90 -o xyz_to_GEOMETRY.x -O3 -static -stand f95
xcart: xyz_to_xcart.f90
	ifort xyz_to_xcart.f90 -o xyz_to_xcart.x -O3 -static -stand f95
clean:
	rm *.x
