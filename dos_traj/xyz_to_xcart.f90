!
!       program raman_traj.f90
!***********************************************************
!
!       
!***********************************************************
!
!

implicit none

integer natom, nsteps, step
integer i, j, k, counter
integer maxline, maxnatom
integer skip, stride
parameter(maxline=100000000)      ! sets a max number of lines that can be read in
parameter(maxnatom=1024)      ! sets a max number of atoms
real(4) snapshot(maxnatom,3)
character(128) fin
character(4) dummyc
character(5) string
character(11) fout

write(6,*) 'Name of data file'
read(5,*) fin

write(6,*) 'Number of configurations to skip?'
read(5,*) skip

write(6,*) 'How large a stride? (1 = all configs used)'
read(5,*) stride

stride = stride - 1

open(1,file=fin,status='old',ERR=100)

read(1,*) natom         ! the first line should be the number of atoms

rewind(1)               ! start at the beginning of the file

counter = 0

do i=1,maxline 
  read(1,*,END=100)
  counter = counter + 1
end do

100 continue

if (counter == maxline)  print *, "Error, did not read to end of file" 

nsteps = counter/(natom + 2)

rewind(1)               ! this puts the read head back at top of file

! skips the initial configurations

counter = 0

do i=1,skip

    read(1,*)
    read(1,*)

    do j=1,natom
        read(1,*)  
    end do

    counter = counter + 1

end do

!

do i=1, (nsteps - skip)/(stride + 1)

    write(unit=string,fmt='(I5.5)') counter + 1 
    fout = "xcart." // string
    open(2, file=fout)

    read(1,*)   
    read(1,*)   

    do j=1,natom
        read(1,*)  dummyc, (snapshot(j,k),k=1,3)
        do k=1,3
            snapshot(j,k) = snapshot(j,k)/0.529177
        end do
        write(2,*) (snapshot(j,k),k=1,3)
    end do

    close(2)

    counter = counter + 1

    do k=1,stride
      read(1,*)
      read(1,*)
      do j=1,natom
          read(1,*)
      end do

      counter = counter + 1

    end do

    if (counter == nsteps) exit

end do 

END
