  SUBROUTINE driver()
    USE io_global,        ONLY : stdout, ionode, ionode_id
    USE io_files,         ONLY : srvaddress
    USE mp_global,        ONLY : mp_startup, mp_bcast, mp_global_end, intra_image_comm
    USE control_flags,    ONLY : conv_elec
    
    USE ions_base,              ONLY : tau
    USE cell_base,              ONLY : alat, at, omega
    USE force_mod,              ONLY : force
    USE ener,                   ONLY : etot
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: MSGLEN=12
    LOGICAL :: isinit=.true., hasdata=.false.
    CHARACTER*12 :: header
    CHARACTER*1024 :: parbuffer
    INTEGER socket, nat
    INTEGER inet, port, ccmd, i
    CHARACTER*1024 :: host  
    REAL *8 :: cellh(3,3), cellih(3,3), vir(3,3), pot
    REAL*8, ALLOCATABLE :: combuf(:)
    REAL*8 :: sigma(3,3)
    inet=1
    host=srvaddress(1:INDEX(srvaddress,':',back=.true.)-1)//achar(0)
    read(srvaddress(INDEX(srvaddress,':',back=.true.)+1:),*) port
    
    IF (ionode) write(*,*) " @ DRIVER MODE: Connecting to host:port ", trim(host), port
    
    IF (srvaddress(1:INDEX(srvaddress,':')-1).eq.('UNIX')) THEN
      inet=0
      host=srvaddress(6:INDEX(srvaddress,':',back=.true.)-1)//achar(0)    
    ENDIF

    IF (ionode) call open_socket(socket, inet, port, host)          
  
    driver_loop: DO
      ! do communication on master node only...
      if (ionode) call readbuffer(socket, header, MSGLEN)      
      call mp_bcast(header,ionode_id, intra_image_comm)
      
      if (ionode) write(*,*) " @ DRIVER MODE: Message from server: ", header
      if (trim(header) == "STATUS") then
         if (ionode) then  ! does not  need init (well, maybe it should, just to check atom numbers and the like... )
            if (hasdata) then
               call writebuffer(socket,"HAVEDATA    ",MSGLEN)
            else
               call writebuffer(socket,"READY       ",MSGLEN)
            endif
         endif
      else if (trim(header) == "POSDATA") then              
         if (ionode) then        
            call readbuffer(socket, cellh, 9*8)
            call readbuffer(socket, cellih, 9*8)
            call readbuffer(socket, nat, 4)
            cellh=transpose(cellh)
            cellih=transpose(cellih)
         endif
         call mp_bcast(cellh,ionode_id, intra_image_comm)
         call mp_bcast(cellih,ionode_id, intra_image_comm)
         call mp_bcast(nat,ionode_id, intra_image_comm)
         if (.not.allocated(combuf)) allocate(combuf(3*nat))
         if (ionode) call readbuffer(socket, combuf, nat*3*8)
         call mp_bcast(combuf,ionode_id, intra_image_comm)
         
         tau = RESHAPE(combuf, (/ 3 , nat /) )/alat                  
         at = cellh / alat
                  
         if (ionode) write(*,*) " @ DRIVER MODE: Received positions "
         CALL hinit1()
         CALL electrons()
         IF ( .NOT. conv_elec ) THEN
           CALL punch( 'all' )
           CALL stop_run( conv_elec )
         ENDIF         
         CALL forces()
         CALL stress(sigma)
         
         combuf=RESHAPE(force, (/ 3 * nat /) ) * 0.5   ! return force in atomic units
         pot=etot * 0.5   ! return potential in atomic units
         vir=transpose(sigma)*omega*0.5   ! return virial in atomic units and without the volume scaling
                  
         hasdata=.true.
      else if (trim(header)=="GETFORCE") then
         if (ionode) write(*,*) " @ DRIVER MODE: Returning v,forces,stress "
         if (ionode) then      
            call writebuffer(socket,"FORCEREADY  ",MSGLEN)            
            call writebuffer(socket,pot,8)
            call writebuffer(socket,nat,4)            
            call writebuffer(socket,combuf,3*nat*8)
            call writebuffer(socket,vir,9*8)
         endif
         hasdata=.false.
         CALL punch( 'config' )
         CALL save_in_ions()
      endif
    ENDDO driver_loop    
    
  END SUBROUTINE
