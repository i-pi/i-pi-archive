      program main
         implicit none

         character*100 :: input_file
         character*8 :: buffer
         character*100 :: output_file
         character*8 :: element
         real :: x, y, z
         integer, parameter :: nat = 108
         integer, parameter :: nsteps = 8000
         integer, parameter :: maxsteps = 4000
         real, parameter :: mass = 2.016*1822.8885
         real, parameter :: atomictime_to_ps = 41341.373
         real, parameter :: bohr_to_angstrom = 1.8897261
         real, parameter :: timestep = 20.0
         real :: p(3,nat,nsteps)
         real :: C(maxsteps)
         real :: mass2nat = 1.0/(mass**2*nat) *
     1                      (atomictime_to_ps/bohr_to_angstrom)**2
         integer :: ncalcs
         integer i, j, k

         call GETARG(1, buffer)

         write(input_file, '(3A)') "test", trim(buffer), ".mom.xyz"
         write(output_file, '(3A)') "mom_corr", trim(buffer), ".out"

         open(unit=11, file=trim(input_file), action="READ")
         open(unit=12, file=trim(output_file), action="WRITE")

         do i = 1, nsteps
            read(11, *) 
            read(11, *) 
            do j = 1, nat
               read(11, '(A8, E13.5, E13.5, E13.5)') element, x, y, z
               p(1,j,i) = x
               p(2,j,i) = y
               p(3,j,i) = z
            end do
            write(*,'(A5 I4 A4 I4 A5)') 
     1                "step ", i, " of ", nsteps, " read"
         end do

         do i = 1, maxsteps
            C(i) = 0.0
            ncalcs = nsteps-i+1
            do j = 1, ncalcs
               do k = 1, nat
                  C(i) = C(i) + dot_product(p(:,k,j),p(:,k,i+j-1))
               end do
            end do
            C(i) = C(i)*mass2nat/ncalcs
            write(*,'(A5 I4 A4 I4 A11)') 
     1                "step ", i, " of ", maxsteps, " calculated"
         end do
         
         do i = 1, maxsteps
            write(12,'(E13.5, E13.5)') 
     1                timestep*(i-1)/atomictime_to_ps, C(i)
            write(*,'(A5 I4 A4 I4 A8)') 
     1                "step ", i, " of ", maxsteps, " printed"
         end do

      end program
