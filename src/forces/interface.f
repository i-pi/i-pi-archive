      program main
         use system_read
         use sys_vars
      implicit none

      character(len=*), parameter :: filedesc = "./system.xml"
      type(Atom), dimension(:), allocatable :: atoms
      type(Cell_vec) :: cell

      call sys_file_read(filedesc, atoms, cell)


      

      end program
