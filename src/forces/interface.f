      program main
         use system_read
         use sys_vars
         use LJ
      implicit none

      character(len=*), parameter :: filedesc = "./system.xml"
      character(len=*), parameter :: filedesc2 = "./system_out.xml"
      type(Atom), dimension(:), allocatable :: atoms
      type(Cell_vec) :: cell
      double precision pot
      double precision, dimension(:,:), allocatable :: f
      double precision, dimension(3,3) :: vir

      integer i, j

      call sys_file_read(filedesc, atoms, cell)
      allocate(f(3,size(atoms)))

      call get_all(atoms, cell, pot, f, vir)
      write(*,'(A, D25.15)') "pot = ", pot 
      write(*,*) "virial:"
      write(*,'(3D25.15)') ((vir(i,j), j = 1, 3), i = 1, 3)
      write(*,*) "forces:"
      write(*,'(3D25.15)') ((f(i,j), i = 1, 3), j = 1, size(atoms))

      call sys_file_write(filedesc2, size(atoms), pot, f, vir)

      deallocate(f)
      deallocate(atoms)

      end program
