      program main
         use system_read
         use sys_vars
         use LJ
      implicit none

      character(len=*), parameter :: filedesc = "./pipepos"
      character(len=*), parameter :: filedesc2 = "./pipeforce"
      type(Atom), dimension(:), allocatable :: atoms
      type(Cell_vec) :: cell
      double precision pot
      double precision, dimension(:,:), allocatable :: f
      double precision, dimension(3,3) :: vir

      integer i, j, ios


      do while (.true.)
         open(555, file = filedesc, iostat = ios, action="read")
         if (ios /= 0) then
            write(*,*) "Error in file reading"
            stop
         end if
         call sys_file_read(555, atoms, cell)
         close(555)
         allocate(f(3,size(atoms)))

         call get_all(atoms, cell, pot, f, vir)

         open(666, file=filedesc2, iostat=ios, action="write")
         if (ios /= 0) then
            write(*,*) "Error in file reading"
            stop
         end if
         call sys_file_write(666, size(atoms), pot, f, vir)
         close(666)
         
         deallocate(f)
         deallocate(atoms)
      enddo

      end program
