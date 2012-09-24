      program main
         implicit none

         character*100 :: input_file
         character*128 :: formatstr
         character*100 :: output_file
         real :: time, corr
         integer, parameter :: maxsteps = 4000
         real, parameter :: timestep = 20.0
         real, parameter :: atomictime_to_ps = 41314.373
         real :: C(maxsteps)
         integer :: ndigits
         integer, parameter :: nruns = 5
         integer, parameter :: ncalcs = 100
         integer :: openfile, nfiles = nruns*ncalcs
         integer i, j, k

         C = 0.0

         write(output_file, '(A16)') "vel_corr_tot.out"
         open(unit=11, file=trim(output_file), action="WRITE")
         do i = 1, nruns
            do j = 1, ncalcs
               ndigits = 1
               do while (j/(10**ndigits) /= 0)
                  ndigits = ndigits + 1
               end do

               write(formatstr, '(A11 I1 A4)') 
     1                  "(A4 I1 A9 I", ndigits, " A4)"
               write(input_file, formatstr) 
     1         "run_", i, "/vel_corr", j, ".out"

               openfile = 11+j+(i-1)*ncalcs
               open(unit=openfile, 
     1                 file=trim(input_file), action="READ")

               do k = 1, maxsteps
                  read(openfile, '(E13.5, E13.5)') time, corr
                  C(k) = C(k) + corr
               end do
               close(openfile)
               write(*,*) 
     1             "read file from run ", i, " and calculation ", j
            end do
         end do

         do k = 1, maxsteps
            write(11,*) timestep*(k-1)/atomictime_to_ps, C(k)/nfiles
            write(*,'(A5 I4 A4 I4 A8)') 
     1                "step ", k, " of ", maxsteps, " printed"
         end do

      end program
