      program main
         use sys_vars
         use LJ
      implicit none

      type(Atom), dimension(:), allocatable :: atoms
      !type(Atom), dimension(:), allocatable :: ref_atoms
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
            ! no need to read stuff, this is not LJ 
            !read(parbuffer(1:nat),*) eps, sigma, rc, rn
            eps = 1.1663e-4
            sigma = 5.270446
            !correction = 4*eps*((sigma/rc)**12 - (sigma/rc)**6)
            correction=0
            rn=rc*1.2
            isinit=.true.
!            write(*,*) "LJ potential initialised with values eps=", 
!     c          eps, ", sigma=",sigma, ", rc=", rc, ", rn=", rn
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

            if ((allocated(n_list)) .neqv. .true.) then
               allocate(f(3,size(atoms)))
               allocate(n_list(size(atoms)-1))
               allocate(index_list((size(atoms)*(size(atoms)-1))/2))
               call nearest_neighbours(atoms, cell, n_list, index_list)
               !allocate(ref_atoms(size(atoms)))
               !ref_atoms = atoms
               allocate(ref_pos(3,size(atoms)))
               do i = 1, nat
                  ref_pos(:,i) = nearest_image(cell, atoms(i)%pos)
               end do
            end if

            do i = 1, size(atoms)
               !q_diff = ref_atoms(i)%pos - atoms(i)%pos
               q_diff = nearest_image(cell,
     1                 nearest_image(cell,atoms(i)%pos) - ref_pos(:,i))
               if (2.0*abs(dot_product(q_diff, q_diff)) >= rn-rc) then
                  call nearest_neighbours(atoms, cell, 
     c                     n_list, index_list)
                  !ref_atoms = atoms
                  !n_list = 0
                  !index_list = 0
                  !counter = 1
                  !do j = 1, size(atoms)-1
                  !   do ios = j+1, size(atoms)
                  !      n_list(j) = n_list(j) + 1
                  !      index_list(counter) = ios
                  !      counter = counter + 1
                  !   end do
                  !end do

                  do j = 1, nat
                     ref_pos(:,j) = nearest_image(cell, atoms(i)%pos)
                  end do
                  exit
               end if
            end do

            call get_all(atoms, cell, n_list, index_list, pot, f, vir)
            vir = transpose(vir)
            ! artificially slow donw evaluation
            ! call sleep(5)
            write(*,*) "computed energy is ",pot
            hasdata=.true.            
         else if (trim(header)=="GETFORCE") then
            call writebuffer(socket,"FORCEREADY  ",MSGLEN)            
            call writebuffer(socket,pot,8)
            call writebuffer(socket,nat,4)            
            do i = 1, nat
               buffer(3*(i-1)+1:3*i)=f(:,i)
            enddo
            ! artificially slow donw transfer
            !call sleep(5)
            call writebuffer(socket,buffer,3*nat*8)
            call writebuffer(socket,vir,9*8)
            hasdata=.false.            
            !hang AFTER transfer is complete
            !call sleep(5)
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

      end program
