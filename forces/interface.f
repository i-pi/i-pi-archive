      program main
         use system_read
         use sys_vars
         use LJ
      implicit none

      type(Atom), dimension(:), allocatable :: atoms
      type(Atom), dimension(:), allocatable :: ref_atoms
      type(Cell_vec) :: cell
      double precision pot
      double precision, dimension(:), allocatable :: buffer
      double precision, dimension(:,:), allocatable :: f
      double precision, dimension(3,3) :: vir
           
      !currently hardcoded, should be 0.2*rc
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
      
      counter = 0
      
      call open_socket(socket)
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
            
            read(parbuffer,*) eps, sigma, rc, rn            
            correction = 4*eps*((sigma/rc)**12 - (sigma/rc)**6)
            isinit=.true.
            write(*,*) "LJ potential initialised with values eps=", 
     c          eps, ", sigma=",sigma, ", rc=", rc, ", rn=", rn
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
               allocate(index_list(size(atoms)*(size(atoms)-1)/2))
               call nearest_neighbours(atoms, cell, n_list, index_list)
               allocate(ref_atoms(size(atoms)))
               ref_atoms = atoms
            end if
            
            do i = 1, size(atoms)
               q_diff = atoms(i)%pos - ref_atoms(i)%pos
               if (2.0*abs(dot_product(q_diff, q_diff)) >= rn-rc) then
                  call nearest_neighbours(atoms, cell, 
     c                     n_list, index_list)
                  ref_atoms = atoms
                  exit
               end if
            end do

            call get_all(atoms, cell, n_list, index_list, pot, f, vir)
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
            hasdata=.false.            
         else
            write(*,*) "Now got ", header
         end if
      enddo
      deallocate(n_list, atoms, f)

      end program
