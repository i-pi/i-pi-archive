      program main
         use system_read
         use sys_vars
         use LJ
      implicit none

      character(len=*), parameter :: filedesc = "../src/pipepos"
      character(len=*), parameter :: filedesc2 = "../src/pipeforce"
      type(Atom), dimension(:), allocatable :: atoms
      type(Atom), dimension(:), allocatable :: ref_atoms
      type(Cell_vec) :: cell
      double precision pot
      double precision, dimension(:,:), allocatable :: f
      double precision, dimension(3,3) :: vir

      double precision, parameter :: skin = 9.33d0
      double precision, dimension(3) :: q_diff

      integer i, j, ios, counter
      logical, dimension(:,:), allocatable :: n_list
      
      counter = 0
      do while (.true.)
         open(555, file = filedesc, iostat = ios, action="read")
         if (ios /= 0) then
            write(*,*) "Error in file reading"
            stop
         end if

         call sys_file_read(555, atoms, cell)
         close(555)

         allocate(f(3,size(atoms)))
!         if ((allocated(n_list)) .neqv. .true.) then
!            allocate(n_list(size(atoms), size(atoms)))
!            call nearest_neighbours(atoms, cell, n_list)
!            allocate(ref_atoms(size(atoms)))
!            ref_atoms = atoms
!         end if
         if ((allocated(n_list)) .neqv. .true.) then
            allocate(n_list(size(atoms), size(atoms)))
         end if
         n_list = .true.

!         do i = 1, size(atoms)
!            q_diff = atoms(i)%pos - ref_atoms(i)%pos
!            if (2.0*abs(dot_product(q_diff, q_diff)) >= skin) then
!               write(*,*) "Calling nearest_neighbours..."
!               call nearest_neighbours(atoms, cell, n_list)
!               ref_atoms = atoms
!               exit
!            end if
!         end do
         call get_all(atoms, cell, n_list, pot, f, vir)

         open(666, file=filedesc2, iostat=ios, action="write")
         if (ios /= 0) then
            write(*,*) "Error in file reading"
            stop
         end if
         call sys_file_write(666, size(atoms), pot, f, vir)
         close(666)
         
         counter = counter + 1
         deallocate(f)
         deallocate(atoms)
      enddo
      deallocate(n_list)

      end program
