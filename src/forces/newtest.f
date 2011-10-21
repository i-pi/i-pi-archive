      program main
         use system_read
         use sys_vars
      implicit none

      character(len=*), parameter :: filedesc = "./system.xml"
      type(Atom), dimension(:), allocatable :: atoms
      type(Cell_vec) :: cell

      integer i, j, k

      call sys_file_read(filedesc, atoms, cell)
      do i = 1, size(atoms)
         write(*,*) (atoms(i)%pos(j), j = 1, 3)
      end do

      write(*,*)

      do j = 1, 3
         write(*,*) (cell%h(i, j), i = 1, 3)
      end do
      do j = 1, 3
         write(*,*) (cell%ih(i, j), i = 1, 3)
      end do

      end program
