      module system_read
         use xml_read
         use sys_vars
      implicit none

      contains

         subroutine sys_file_read(filedesc, atoms, cell)
            character(len=*), intent(in) :: filedesc
            type(Atom), dimension(:), allocatable, intent(out) 
     1:: atoms
            type(Cell_vec), intent(out) :: cell

            character(len=200) :: file_line
            integer natoms, ios, counter, i, j, k
            logical correct
            double precision, dimension(3) :: temp_array

            if (allocated(atoms)) deallocate(atoms)

            open(5, file = filedesc, iostat = ios)
            if (ios /= 0) then
               write(*,*) "Error in file reading"
               stop
            end if

            read(5,'(A200)')
            read(5,'(A200)') file_line
            counter = 2
            call search_begin(file_line, "System", correct)
            if (.not. correct) then
               write(*,*) "Error in line 2, rootname not System"
               stop
            end if

            read(5,'(A200)') file_line
            counter = counter + 1
            call read_value(file_line, "natoms", natoms, correct)
            if (.not. correct) then
               write(*,*) "Error in line 3"
               stop
            else
               allocate(atoms(natoms))
            end if

            do i = 1, natoms
               read(5,'(A200)') file_line
               counter = counter + 1
               call search_begin(file_line, "Atom_vec", correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  stop
               end if

               read(5,'(A200)') file_line
               counter = counter + 1
               call read_value(file_line, "q", temp_array, correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  stop
               else
                  atoms(i)%pos = temp_array
               end if

               read(5,'(A200)') file_line
               counter = counter + 1
               call search_end(file_line, "Atom_vec", correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  stop
               end if

            end do

            do i = 1, 3
               read(5,'(A200)') file_line
               counter = counter + 1
               call search_begin(file_line, "Cell_vec", correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  stop
               end if

               read(5,'(A200)') file_line
               counter = counter + 1
               call read_value(file_line, "h", temp_array, correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  stop
               else
                  cell%h(:,i) = temp_array
               end if

               read(5,'(A200)') file_line
               counter = counter + 1
               call search_end(file_line, "Cell_vec", correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  stop
               end if
            end do

            do i = 1, 3
               read(5,'(A200)') file_line
               counter = counter + 1
               call search_begin(file_line, "Cell_vec", correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  stop
               end if

               read(5,'(A200)') file_line
               counter = counter + 1
               call read_value(file_line, "ih", temp_array, correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  stop
               else
                  cell%ih(:,i) = temp_array
               end if

               read(5,'(A200)') file_line
               counter = counter + 1
               call search_end(file_line, "Cell_vec", correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  stop
               end if
            end do

            read(5,'(A200)') file_line
            counter = counter + 1
            call search_end(file_line, "System", correct)
            if (.not. correct) then
               write(*,'(A14, I4)') "Error in line ", counter
               stop
            end if

         end subroutine
      
      end module
