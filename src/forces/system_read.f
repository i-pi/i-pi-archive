      module system_read
         use xml_read
         use sys_vars
      implicit none

      contains

         subroutine sys_file_read(filedesc, atoms, cell)
            integer, intent(in) :: filedesc
            type(Atom), dimension(:), allocatable, intent(out) 
     1:: atoms
            type(Cell_vec), intent(out) :: cell

            character(len=200) :: file_line
            integer natoms, counter, i, j, k
            logical correct
            double precision, dimension(3) :: temp_array

            if (allocated(atoms)) deallocate(atoms)

            read(filedesc,'(A200)')
            read(filedesc,'(A200)') file_line
            counter = 2
            write(*,*) "here we are 1"
            call search_begin(file_line, "System", correct)
            if (.not. correct) then
               write(*,*) "Error in line 2, rootname not System"
               stop
            end if
            write(*,*) "here we are 2"
            read(filedesc,'(A200)') file_line
            counter = counter + 1
            call read_value(file_line, "natoms", natoms, correct)
            if (.not. correct) then
               write(*,*) "Error in line 3"
               stop
            else
               allocate(atoms(natoms))
            end if
            write(*,*) "here we are 3"
            do i = 1, natoms
               write(*,*) "Reading atom", i
               read(filedesc,'(A200)') file_line
               counter = counter + 1
               call search_begin(file_line, "Atom_vec", correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  deallocate(atoms)
                  stop
               end if

               read(filedesc,'(A200)') file_line
               counter = counter + 1
               call read_value(file_line, "q", temp_array, correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  deallocate(atoms)
                  stop
               else
                  atoms(i)%pos = temp_array
               end if

               read(filedesc,'(A200)') file_line
               counter = counter + 1
               call search_end(file_line, "Atom_vec", correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  deallocate(atoms)
                  stop
               end if
            end do
            write(*,*) "here we are 4 -- good to go"
            do i = 1, 3
               read(filedesc,'(A200)') file_line
               counter = counter + 1
               call search_begin(file_line, "Cell_vec", correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  deallocate(atoms)
                  stop
               end if

               read(filedesc,'(A200)') file_line
               counter = counter + 1
               call read_value(file_line, "h", temp_array, correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  deallocate(atoms)
                  stop
               else
                  cell%h(:,i) = temp_array
               end if

               read(filedesc,'(A200)') file_line
               counter = counter + 1
               call search_end(file_line, "Cell_vec", correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  deallocate(atoms)
                  stop
               end if
            end do
            write(*,*) "here we are 5 -- good to go"
            do i = 1, 3
               write(*,*) "reading line ", i
               read(filedesc,'(A200)') file_line
               counter = counter + 1
               call search_begin(file_line, "Cell_vec", correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  deallocate(atoms)
                  stop
               end if
               
               write(*,*) "reading HERE ", i
               read(filedesc,'(A200)') file_line               
               counter = counter + 1
               call read_value(file_line, "ih", temp_array, correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  deallocate(atoms)
                  stop
               else
                  cell%ih(:,i) = temp_array
               end if

              write(*,*) "reading HERE2 ", file_line
 
               read(filedesc,'(A)') file_line
              write(*,*) "read HERE2 ", file_line

               counter = counter + 1
               call search_end(file_line, "Cell_vec", correct)
               if (.not. correct) then
                  write(*,'(A14, I4)') "Error in line ", counter
                  deallocate(atoms)
                  stop
               end if
            end do
            write(*,*) "here we are 6 -- good to go"

            read(filedesc,'(A200)') file_line
            counter = counter + 1
            call search_end(file_line, "System", correct)
            if (.not. correct) then
               write(*,'(A14, I4)') "Error in line ", counter
               deallocate(atoms)
               stop
            end if
            write(*,*) "here we are -- really done"
         end subroutine

         subroutine sys_file_write(filedesc, natoms, pot, f, vir)
            integer, intent(in) :: filedesc
            double precision, intent(in) :: pot
            integer, intent(in) :: natoms
            double precision, dimension(3,natoms), intent(in) :: f
            double precision, dimension(3,3), intent(in) :: vir

            integer ios
            character(len=200) :: file_line
            character(len=3), parameter :: tab = "    "

            integer i


            write(filedesc,'(A)') "<?xml version='1.0'?>"

            call write_begin("System", file_line)    
            write(filedesc,'(A)') trim(file_line)

            call write_value("pot", file_line, pot)
            write(filedesc, '(A3, A)') tab, trim(file_line)

            do i = 1, natoms
               call write_begin("atom_f", file_line)
               write(filedesc, '(A3, A)') tab, trim(file_line) 
               
               call write_value("f", file_line, f(:,i))
               write(filedesc,'(A3, A3, A)') tab, tab, trim(file_line)

               call write_end("atom_f", file_line)
               write(filedesc, '(A3, A)') tab, trim(file_line) 
            end do

            do i = 1, 3
               call write_begin("vir_column", file_line)
               write(filedesc, '(A3, A)') tab, trim(file_line) 
               
               call write_value("x", file_line, vir(:,i))
               write(filedesc,'(A3, A3, A)') tab, tab, trim(file_line)

               call write_end("vir_column", file_line)
               write(filedesc, '(A3, A)') tab, trim(file_line) 
            end do

            call write_end("System", file_line)    
            write(filedesc,'(A)') trim(file_line)

         end subroutine
      
      end module
