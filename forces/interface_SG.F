      program main
         use sys_vars
         use SG
      implicit none

      type(Atom), dimension(:), allocatable :: atoms
      type(Cell_vec) :: cell
      double precision, dimension(:,:), allocatable :: ref_pos
      double precision pot
      double precision, dimension(:), allocatable :: buffer
      double precision, dimension(:,:), allocatable :: f
      double precision, dimension(3,3) :: vir
           
      double precision, dimension(3) :: q_diff

      integer i, j, ios, counter
      integer, dimension(:), allocatable :: n_list
      integer, dimension(:), allocatable :: index_list
      double precision :: time = 0.0, timeall=0.0, timewait=0.0
      integer countnum, countrate
      
      integer, parameter :: MSGLEN=12
      logical :: isinit=.false., hasdata=.false.
      character*12 :: header
      character*1024 :: parbuffer
      integer socket, nat
      integer*4 inet, port, ccmd
      character*1024 :: host
      
      ccmd=0
      inet=1
      host="localhost"//achar(0)
      port=31415
      do i=1, IARGC()
         call GETARG(i, parbuffer)
         if (parbuffer == "-u") then 
            inet=0
            ccmd=0
         elseif (parbuffer == "-h") then 
            ccmd=1
         elseif (parbuffer == "-p") then 
            ccmd=2
         elseif (parbuffer == "-r") then 
            ccmd=3
         else
            if (ccmd==0) then 
               write(*,*) "Unrecognized command line argument", ccmd
               write(*,*) 
     c "SYNTAX: ljengine.x [-u] -h hostname -p port -r cutoff"
               call exit(-1)
            endif
            if (ccmd==1) then
               host=trim(parbuffer)//achar(0)
            elseif (ccmd==2) then
               read(parbuffer,*) port
            elseif (ccmd==3) then
               read(parbuffer,*) rc
            endif
            ccmd=0            
         endif
      enddo
      write(*,*) "Connecting with options  ", trim(host), port, inet
      counter = 0
      
      call open_socket(socket, inet, port, host)
      do while (.true.)
         
         call readbuffer(socket, header, MSGLEN)
         write(*,*) "Message from server: ", header

         if (trim(header) == "STATUS") then
            if (.not. isinit) then
               call writebuffer(socket,"NEEDINIT    ",MSGLEN)
            else if (hasdata) then
               call writebuffer(socket,"HAVEDATA    ",MSGLEN)
            else
               call writebuffer(socket,"READY       ",MSGLEN)
            endif

         else if (trim(header) == "INIT") then     
            call readbuffer(socket, nat, 4)
            call readbuffer(socket, parbuffer, nat)
            correction=0
            rn=rc*1.2
!            correction = exp(1.713d0-1.5671d0*rc-0.00993d0*rc**2) -
!     1 (12.14d0/rc**6 + 215.2d0/rc**8 - 143.1d0/rc**9 + 4813.9d0/rc**10)
            isinit=.true.

         else if (trim(header) == "POSDATA") then              
            call readbuffer(socket, cell%h, 9*8)
            call readbuffer(socket, cell%ih, 9*8)
            cell%h=transpose(cell%h)    !fortran data arrays held transposed
            cell%ih=transpose(cell%ih)  !compared to python arrays
            call readbuffer(socket, nat, 4)
            if ( .not. allocated(buffer) ) then
               write(*,*) "allocating buffer"
               allocate(buffer(3*nat))
            endif
            call readbuffer(socket, buffer, nat*3*8)
            
            if ((allocated(atoms)) .neqv. .true.) then
               write(*,*) "allocating nlist"      
               allocate(atoms(nat))      
            end if
            do i = 1, nat
               atoms(i)%pos=buffer(3*(i-1)+1:3*i)
            enddo

            if ((allocated(n_list)) .neqv. .true.) then
               allocate(f(3,size(atoms)))
               allocate(n_list(size(atoms)-1))
               allocate(index_list((size(atoms)*(size(atoms)-1))/2))
               call nearest_neighbours(atoms, cell, n_list, index_list)
               allocate(ref_pos(3,size(atoms)))
               do i = 1, nat
                  ref_pos(:,i) = atoms(i)%pos
               end do
            end if

            do i = 1, size(atoms)
               q_diff = nearest_image(cell,
     1                 atoms(i)%pos - ref_pos(:,i))

               if (2.0*sqrt(dot_product(q_diff, q_diff)) >= rn-rc) then
                  call nearest_neighbours(atoms, cell, 
     c                     n_list, index_list)

                  do j = 1, nat
                     ref_pos(:,j) = atoms(j)%pos
                  end do
                  exit
               end if
            end do

            call get_all(atoms, cell, n_list, index_list, pot, f, vir)
            call long_range_corr(cell, nat, vir, pot)
            vir = transpose(vir)
            write(*,*) "computed energy is ",pot
            hasdata=.true.            

         else if (trim(header)=="GETFORCE") then
            call writebuffer(socket,"FORCEREADY  ",MSGLEN)            
            call writebuffer(socket,pot,8)
            call writebuffer(socket,nat,4)            
            do i = 1, nat
               buffer(3*(i-1)+1:3*i)=f(:,i)
            enddo
            call writebuffer(socket,buffer,3*nat*8)
            call writebuffer(socket,vir,9*8)
!            nat=12
            nat=0
            call writebuffer(socket,nat,4)
!            call writebuffer(socket,"abcdefghijkl",nat)            
            hasdata=.false.            
         else
            write(*,*) "Now got ", header
         end if

      enddo

      deallocate(n_list, atoms, f)

      contains

         function nearest_image(cell, q)
            Type(Cell_vec), intent(in) :: cell
            double precision, dimension(3), intent(in) :: q

            double precision, dimension(3) :: nearest_image
            double precision, dimension(3) :: s

            integer i

            s = matmul(cell%ih, q)
            do i = 1, 3
               s(i) = s(i) - dnint(s(i))
            end do
            s = matmul(cell%h, s)
            do i = 1, 3
               nearest_image(i) = s(i)
            end do

         end function   

         subroutine long_range_corr(cell, natoms, stress, pot)
            Type(Cell_vec), intent(in) :: cell
            integer, intent(in) :: natoms

            double precision, dimension(3,3), intent(inout) :: stress
            double precision, intent(inout) :: pot

            double precision :: volume, pot_corr, pi, normalization
            double precision :: dpot, dforce

            pi = 3.14159265358979

            volume = cell%h(1,1)*cell%h(2,2)*cell%h(3,3)

            normalization = 2*pi*natoms**2/volume

            pot_corr = normalization*long_range(rc) 

            pot = pot + pot_corr

            call SG_functions(rc, dpot, dforce)
            do i = 1,3
               stress(i,i) = stress(i,i) + (normalization/3.0*rc**3
     1                          *dpot + pot_corr)
            end do

         end subroutine

      end program
