PROGRAM dos_kpt 
!
!       program dos_kpt.f90
!***********************************************************
!
!       
!***********************************************************
!

implicit none

integer i, j, k, nkpt, neig, nedos, diagonal_n, nstates
integer maxline, counter
parameter(maxline=10000000)      ! this just sets a max number of lines that can be read in
real(8) E, emin, emax, dE, weight, eigenvalue, sigma
real(8) dft_sum 
real(8), allocatable :: DFT_dos(:,:)
real(8), allocatable :: wtk(:)
real(8), allocatable :: DFT_eig_array(:)
real, parameter :: pi = 3.14159 ! cannot be changed
character(128) fin, fin2
character(6) dummya, dummyb, dummyc, dummyd, dummye 

nedos = 1101 
sigma = 0.1

fin = 'band.out'
open(1,file=fin,status='old',ERR=100)

do i=1,maxline 
  read(1,*,END=100)
  counter = counter + 1
end do

100  continue

if (counter == maxline)  print *, "Error, did not read to end of file" 

nstates = 1
diagonal_n = counter - 2 

rewind(1)               ! this puts the read head back at top of file
read(1,*)               ! skips the header

allocate(DFT_dos(3,nedos))
allocate(wtk(nstates*diagonal_n))
allocate(DFT_eig_array(nstates*diagonal_n))

DFT_dos(:,:) = 0.0
DFT_eig_array(:) = 0.0
wtk(:) = 1.0

do j=1,nstates*diagonal_n
       read(1,*) DFT_eig_array(j)
end do

emin = -40
emax = 10 


dE = (emax - emin)/nedos

do i=1,nedos

    E = (i - 1)*dE + emin  ! should I change this to middle..?
                           ! no, this is the left side of the box for
                           ! integration

    do j=1,nstates*diagonal_n

            weight = 1.0*wtk(j)  ! two e per state
            eigenvalue = DFT_eig_array(j)
            DFT_dos(2,i) = DFT_dos(2,i) + weight*(1.0/2.0)*(erf(sqrt(2.0)*(E + dE -eigenvalue)/(2*sigma)) &
            - erf(sqrt(2.0)*(E - eigenvalue)/(2*sigma)))

    end do

end do

open(3,file="dft_dos.dat")

dft_sum = 0.0

do i=1,nedos
  E = (i-1)*dE + emin + dE/2
  dft_sum = dft_sum + DFT_dos(2,i)
  write(3,*) E, DFT_dos(2,i), dft_sum
end do

close(3)

deallocate(wtk)
deallocate(DFT_eig_array)
deallocate(DFT_dos)

END PROGRAM dos_kpt 

