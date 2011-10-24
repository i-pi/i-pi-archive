      module xml_read
      implicit none

      interface read_value
         module procedure read_array, read_real, read_integer
      end interface
      interface write_value
         module procedure write_array, write_real, write_integer
      end interface

      private
      public :: read_value, search_begin, search_end
      public :: write_begin, write_end, write_value

      contains

         subroutine search_begin(file_line, tag, correct)
            character(len=*), intent(in) :: tag
            character(len=*), intent(in) :: file_line
            logical, intent(out) :: correct

            integer first, last
            character(len=len(file_line)) :: input_tag

            first = scan(file_line, "<")
            last = scan(file_line, ">")

            if (first == 0 .or. last == 0 .or. last <= first+1) then
               write(*,*) "Tag header not formatted correctly in ", tag
               correct = .false.
            end if

            input_tag = file_line(first+1:last-1)
            if (trim(input_tag) == tag) then
               correct = .true.
            else
               write(*,*) "Incorrect tag in ", tag
               correct = .false.
            end if

         end subroutine

         subroutine search_end(file_line, tag, correct)
            character(len=*), intent(in) :: tag
            character(len=*), intent(in) :: file_line
            logical, intent(out) :: correct

            integer first, second, last
            character(len=len(file_line)) :: input_tag

            first = scan(file_line, "<", .true.)
            second = scan(file_line, "/", .true.)
            last = scan(file_line, ">", .true.)

            if (first == 0 .or. last == 0 .or. last <= second+1 .or. 
     1second /= first+1) then
               write(*,*) "Tag closer not formatted correctly in ", tag
               correct = .false.
            end if

            input_tag = file_line(second+1:last-1)
            if (trim(input_tag) == tag) then
               correct = .true.
            else
               write(*,*) "Incorrect tag in ", tag
               correct = .false.
            end if

         end subroutine

         logical function format_test(format_vars)
            integer, dimension(6), intent(in) :: format_vars
         
            integer i

            do i = 1,6
               if (format_vars(i) == 0) then
                  write(*,*) "Missing format labels in array"
                  format_test = .false.
                  return
               end if
            end do
            
            if (format_vars(2) /= format_vars(1)+1 .or. 
     1format_vars(6) /= format_vars(5)+1) then
               write(*,*) "Extra symbols between the array and tag"
               format_test = .false.
               return
            end if

            if (format_vars(3) == format_vars(2)+1 .or.
     1format_vars(4) == format_vars(3)+1 .or. format_vars(5) ==
     2format_vars(4)+1) then
               write(*,*) "Missing elements in array"
               format_test = .false.
               return
            end if

            do i = 1,5
               if (format_vars(i) >= format_vars(i+1)) then
                  write(*,*) "Format error in array"
                  format_test = .false.
                  return
               else
                  format_test = .true.
                  return
               end if
            end do

         end function

         subroutine read_array(file_line, tag, array, success)
            character(len=*), intent(in) :: tag
            character(len=*), intent(in) :: file_line
            double precision, dimension(3), intent(out) :: array
            logical, intent(out) :: success

            integer ios1, ios2, ios3
            logical correct_begin, correct_end, correct_total
            integer, dimension(6) :: format_vars

            ios1 = 0
            ios2 = 0
            ios3 = 0
            success = .false.

            call search_begin(file_line, tag, correct_begin)
            call search_end(file_line, tag, correct_end)

            if (correct_begin .and. correct_end) then
               format_vars(1) = scan(file_line, ">")
               format_vars(2) = scan(file_line, "[")
               format_vars(3) = scan(file_line, ",")
               format_vars(4) = scan(file_line, ",", .true.)
               format_vars(5) = scan(file_line, "]", .true.)
               format_vars(6) = scan(file_line, "<", .true.)

               correct_total = format_test(format_vars)
               
               if (correct_total) then
                  read(file_line(format_vars(2)+1:format_vars(3)-1),
     1'(D20.15)', iostat = ios1) array(1)
                  read(file_line(format_vars(3)+1:format_vars(4)-1),
     1'(D20.15)', iostat = ios2) array(2)
                  read(file_line(format_vars(4)+1:format_vars(5)-1),
     1'(D20.15)', iostat = ios3) array(3)

                  if (ios1 /= 0 .or. ios2 /= 0 .or. ios3 /= 0) then
                     write(*,*) "Tried to write NaN to array in ", tag
                     array = 0.0
                  else
                     success = .true.
                  end if
               else
                  write(*,*) "Error in tag ", tag
                  array = 0.0
               end if
            else
               write(*,*) "Error in tag ", tag
               array = 0.0
            end if

         end subroutine

         subroutine read_real(file_line, tag, value, success)
            character(len=*), intent(in) :: tag
            character(len=*), intent(in) :: file_line
            double precision, intent(out) :: value
            logical, intent(out) :: success

            integer first, last, ios
            logical correct_begin, correct_end

            ios = 0
            success = .false.

            call search_begin(file_line, tag, correct_begin)
            call search_end(file_line, tag, correct_end)

            if (correct_begin .and. correct_end) then
               first = scan(file_line, ">")
               last = scan(file_line, "<", .true.)

               if (first > last) then
                  write(*,*) "Format error in tag ", tag
                  value = 0.0
               else if (first == last-1) then
                  write(*,*) "Missing number between tags in ", tag
                  value = 0.0
               else
                  read(file_line(first+1:last-1), '(D20.15)', 
     1iostat = ios) value
                  if (ios /= 0) then
                     write(*,*) "Tried to write NaN to real in ", tag
                     value = 0.0
                  else
                     success = .true.
                  end if
               end if
            else
               write(*,*) "Error in tag ", tag
               value = 0.0
            end if

         end subroutine

         subroutine read_integer(file_line, tag, value, success)
            character(len=*), intent(in) :: tag
            character(len=*), intent(in) :: file_line
            integer, intent(out) :: value
            logical, intent(out) :: success

            integer first, last, ios
            logical correct_begin, correct_end

            ios = 0
            success = .false.

            call search_begin(file_line, tag, correct_begin)
            call search_end(file_line, tag, correct_end)

            if (correct_begin .and. correct_end) then
               first = scan(file_line, ">")
               last = scan(file_line, "<", .true.)

               if (first > last) then
                  write(*,*) "Format error in tag ", tag
                  value = 0
               else if (first == last-1) then
                  write(*,*) "Missing number between tags in ", tag
                  value = 0
               else
                  read(file_line(first+1:last-1), '(I10)', 
     1iostat = ios) value
                  if (ios /= 0) then
                     write(*,*) "Tried to write NaN to integer in ", tag
                     value = 0
                  else
                     success = .true.
                  end if
               end if
            else
               write(*,*) "Error in tag ", tag
               value = 0
            end if

         end subroutine

         subroutine write_begin(tag, file_line)
            character(len=*), intent(in) :: tag
            character(len=200), intent(out) :: file_line

            write(file_line,*) "<", tag, ">"

         end subroutine

         subroutine write_end(tag, file_line)
            character(len=*), intent(in) :: tag
            character(len=200) :: file_line

            write(file_line,*) "</", tag, ">"

         end subroutine

         subroutine write_real(tag, file_line, value)
            character(len=*), intent(in) :: tag
            character(len=200), intent(out) :: file_line
            double precision, intent(in) :: value

            write(file_line,'(3A, D25.15, 3A)') 
     1"<", tag, ">", value, "</", tag, ">"
   
         end subroutine

         subroutine write_integer(tag, file_line, value)
            character(len=*), intent(in) :: tag
            character(len=200), intent(out) :: file_line
            integer, intent(in) :: value

            write(file_line,'(3A, I10, 3A)') 
     1"<", tag, ">", value, "</", tag, ">"
   
         end subroutine

         subroutine write_array(tag, file_line, value)
            character(len=*), intent(in) :: tag
            character(len=200), intent(out) :: file_line
            double precision, dimension(3), intent(in) :: value

            write(file_line, '(3A, D25.15, A1, D25.15, A1, D25.15, 3A)')
     1"<", tag, ">[", value(1), ",", value(2), ",", value(3),"]</",
     2tag, ">"

         end subroutine

      end module
