      program main

      use sys_vars
      implicit none

      type(Atom), dimension(:), allocatable :: atoms
      type(Cell_vec) :: cell
      double precision, dimension(:,:), allocatable :: ref_pos
      double precision pot, kq
      double precision, dimension(:), allocatable :: buffer
      double precision, dimension(:,:), allocatable :: f
      double precision, dimension(3,3) :: vir

      double precision, dimension(3) :: dx
      double precision d, du

      integer i, j, k, ios, counter
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
         elseif (parbuffer == "-k") then
            ccmd=3
         else
            if (ccmd==0) then
               write(*,*) "Unrecognized command line argument", ccmd
               write(*,*)
     c "SYNTAX: gas_mol.x [-u] -h hostname -p port -r cutoff -k spring"
               call exit(-1)
            endif
            if (ccmd==1) then
               host=trim(parbuffer)//achar(0)
            elseif (ccmd==2) then
               read(parbuffer,*) port
            elseif (ccmd==3) then
               read(parbuffer,*) kq
            endif
            ccmd=0
         endif
      enddo
      write(*,*) "Connecting with options  ", trim(host), port, inet
      counter = 0

      call open_socket(socket, inet, port, host)
      do while (.true.)

         call readbuffer(socket, header, MSGLEN)
         !write(*,*) "Message from server: ", header
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
            isinit=.true.
         else if (trim(header) == "POSDATA") then
            call readbuffer(socket, cell%h, 9*8)
            call readbuffer(socket, cell%ih, 9*8)
            cell%h=transpose(cell%h)
            cell%ih=transpose(cell%ih)
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

            if ((allocated(f)) .neqv. .true.) then
               allocate(f(3,size(atoms)))
            end if


            dx = atoms(2)%pos-atoms(1)%pos
            d = dot_product(dx,dx)
            pot=0.5*kq*d
            d = sqrt(d)
            du = kq*d

            f(:,1)= du*dx/d
            f(:,2)= -du*dx/d


            vir = 0.0
            do i = 1,nat
               do j = 1,3
                  do k=1,3
                   vir(j,k) = vir(j,k) + f(j,i)*dx(k)
                  enddo
               enddo
            enddo

            vir = transpose(vir)

            ! artificially slow down evaluation
            !call sleep(2)
            !write(*,*) "computed energy is ",pot
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
            nat=0
            call writebuffer(socket,nat,4)
            hasdata=.false.
         else
            write(*,*) "Now got ", header
         end if
      enddo
      deallocate( atoms, f)

      end program
