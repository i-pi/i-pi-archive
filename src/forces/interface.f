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
         write(*,*) "opened reading pipe"
         if (ios /= 0) then
            write(*,*) "Error in file reading"
            stop
         end if
         call sys_file_read(555, atoms, cell)
         write(*,*) "file has been read"
         close(555)
         allocate(f(3,size(atoms)))

         call get_all(atoms, cell, pot, f, vir)
         write(*,'(A, D25.15)') "pot = ", pot 
         write(*,*) "virial:"
         write(*,'(3D25.15)') ((vir(i,j), j = 1, 3), i = 1, 3)
         write(*,*) "forces:"
         write(*,'(3D25.15)') ((f(i,j), i = 1, 3), j = 1, size(atoms))


         open(666, file=filedesc2, iostat=ios, action="write")
         write(*,*) "opened writing pipe"
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
