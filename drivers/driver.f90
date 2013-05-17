      PROGRAM DRIVER
         USE LJ
         USE SG
      IMPLICIT NONE

      ! SOCKET VARIABLES
      INTEGER, PARAMETER :: MSGLEN=12   ! length of the headers of the driver/wrapper communication protocol
      INTEGER socket, inet, port        ! socket ID & address of the server
      CHARACTER*1024 :: host

      ! COMMAND LINE PARSING
      CHARACTER*1024 :: cmdbuffer, vops
      INTEGER ccmd, vstyle
      LOGICAL verbose
      INTEGER commas(4), par_count      ! stores the index of commas in the parameter string
      DOUBLE PRECISION vpars(4)         ! array to store the parameters of the potential

      ! SOCKET COMMUNICATION BUFFERS
      CHARACTER*12 :: header
      LOGICAL :: isinit=.false., hasdata=.false.
      INTEGER cbuf
      CHARACTER*2048 :: initbuffer      ! it's unlikely a string this large will ever be passed...
      DOUBLE PRECISION, ALLOCATABLE :: msgbuffer(:)

      ! PARAMETERS OF THE SYSTEM (CELL, ATOM POSITIONS, ...)
      INTEGER nat
      DOUBLE PRECISION pot
      DOUBLE PRECISION, ALLOCATABLE :: atoms(:,:), forces(:,:)
      DOUBLE PRECISION cell_h(3,3), cell_ih(3,3), virial(3,3)


      integer i, j, k, ios, counter
      integer, dimension(:), allocatable :: n_list
      integer, dimension(:), allocatable :: index_list
      double precision :: time = 0.0, timeall=0.0, timewait=0.0
      integer countnum, countrate


      ! parse the command line parameters
      ! intialize defaults
      ccmd = 0
      inet = 1
      host = "localhost"//achar(0)
      port = 31415
      verbose = .false.
      par_count = 0

      DO i=1, IARGC()
         CALL GETARG(i, cmdbuffer)
         IF (cmdbuffer == "-u") THEN
            inet=0
            ccmd=0
         ELSEIF (cmdbuffer == "-h") THEN
            ccmd=1
         ELSEIF (cmdbuffer == "-p") THEN
            ccmd=2
         ELSEIF (cmdbuffer == "-m") THEN ! reads the style of the potential function
            ccmd=3
         ELSEIF (cmdbuffer == "-o") THEN
            ccmd=4
         ELSEIF (cmdbuffer == "-v") THEN
            verbose = .true.
         ELSE
            IF (ccmd==0) THEN
               WRITE(*,*) " Unrecognized command line argument", ccmd
               WRITE(*,*) " SYNTAX: driver.x [-u] -h hostname -p port -m [gas|lj|sg] -o 'comma_separated_parameters' [-v] "
               WRITE(*,*) ""
               WRITE(*,*) " For LJ potential use -o sigma,epsilon,cutoff "
               WRITE(*,*) " For SG potential use -o cutoff "
               WRITE(*,*) " For the ideal gas, no options needed! "
               CALL EXIT(-1)
            ENDIF
            IF (ccmd == 1) THEN
               host = trim(cmdbuffer)//achar(0)
            ELSEIF (ccmd == 2) THEN
               READ(cmdbuffer,*) port
            ELSEIF (ccmd == 3) THEN
               IF (trim(cmdbuffer) == "lj") THEN
                  vstyle = 1
               ELSEIF (trim(cmdbuffer) == "sg") THEN
                  vstyle = 2
               ELSE
                  vstyle = 0  ! ideal gas
               ENDIF
            ELSEIF (ccmd == 4) THEN
               par_count = 1
               commas(1) = 0
               DO WHILE (index(cmdbuffer(commas(par_count)+1:), ',') > 0) 
                  commas(par_count + 1) = index(cmdbuffer(commas(par_count)+1:), ',') + commas(par_count)
                  READ(cmdbuffer(commas(par_count)+1:commas(par_count + 1)-1),*) vpars(par_count)
                  par_count = par_count + 1
               ENDDO
               READ(cmdbuffer(commas(par_count)+1:),*) vpars(par_count)
            ENDIF
            ccmd=0
         ENDIF
      ENDDO

      IF (verbose) THEN
         WRITE(*,*) " DRIVER - Connecting to host ", trim(host)
         IF (inet > 0) THEN
            WRITE(*,*) " on port ", port, " using an internet socket."
         ELSE
            WRITE(*,*) " using an UNIX socket."
         ENDIF
      ENDIF

      ! Calls the interface to the C sockets to open a communication channel
      CALL open_socket(socket, inet, port, host)
      nat = -1
      DO WHILE (.true.) ! Loops forever (or until the wrapper ends!)

         ! Reads from the socket one message header
         CALL readbuffer(socket, header, MSGLEN)
         IF (verbose) WRITE(*,*) " Message from server: ", trim(header)

         IF (trim(header) == "STATUS") THEN
            ! The wrapper is inquiring on what we are doing
            IF (.not. isinit) THEN
               CALL writebuffer(socket,"NEEDINIT    ",MSGLEN)  ! Signals that we need initialization data
            ELSEIF (hasdata) THEN
               CALL writebuffer(socket,"HAVEDATA    ",MSGLEN)  ! Signals that we are done computing and can return forces
            ELSE
               CALL writebuffer(socket,"READY       ",MSGLEN)  ! We are idling and eager to compute something
            ENDIF
         ELSEIF (trim(header) == "INIT") THEN     ! The driver is kindly providing a string for initialization
            CALL readbuffer(socket, cbuf, 4)
            CALL readbuffer(socket, initbuffer, cbuf)
            IF (verbose) WRITE(*,*) " Initializing system from wrapper, using ", trim(initbuffer)
            isinit=.true. ! We actually do nothing with this string, thanks anyway. Could be used to pass some information (e.g. the input parameters, or the index of the replica, from the driver
         ELSEIF (trim(header) == "POSDATA") THEN  ! The driver is sending the positions of the atoms. Here is where we do the calculation!

            ! Parses the flow of data from the socket
            CALL readbuffer(socket, cell_h,  9*8)  ! Cell matrix
            CALL readbuffer(socket, cell_ih, 9*8)  ! Inverse of the cell matrix (so we don't have to invert it every time here)
            CALL readbuffer(socket, cbuf, 4)       ! The number of atoms in the cell

            IF (nat < 0) THEN  ! Assumes that the number of atoms does not change throughout a simulation
               nat = cbuf
               IF (verbose) WRITE(*,*) " Allocating buffers, with ", nat, " atoms"
               ALLOCATE(msgbuffer(3*nat))
               ALLOCATE(atoms(3,nat))
               ALLOCATE(forces(3,nat))
            ENDIF

            CALL readbuffer(socket, msgbuffer, nat*3*8)

            ! The wrapper uses atomic units for everything, and row major storage.
            ! At this stage one should take care that everything is converted in the
            ! units and storage mode used in the driver.
            cell_h = transpose(cell_h)
            cell_ih = transpose(cell_ih)

            atoms = reshape(msgbuffer,  (/ 3, nat /) )

            IF (vstyle == 0) THEN   ! ideal gas, so no calculation done
               pot = 0
               forces = 0
               virial = 0
            ENDIF
            hasdata=.true. ! Signal that we have data ready to be passed back to the wrapper
         ELSEIF (trim(header) == "GETFORCE") THEN  ! The driver calculation is finished, it's time to send the results back to the wrapper

            ! Data must be re-formatted (and units converted) in the units and shapes used in the wrapper
            msgbuffer = reshape(forces, (/ 3*nat /) )
            virial = transpose(virial)

            CALL writebuffer(socket,"FORCEREADY  ",MSGLEN)
            CALL writebuffer(socket,pot,8)  ! Writing the potential
            CALL writebuffer(socket,nat,4)  ! Writing the number of atoms
            CALL writebuffer(socket,msgbuffer,3*nat*8) ! Writing the forces
            CALL writebuffer(socket,virial,9*8)  ! Writing the virial tensor, NOT divided by the volume
            cbuf = 7
            CALL writebuffer(socket,cbuf,4) ! This would write out the "extras" string, but in this case there isn't one.
            CALL writebuffer(socket,"nothing",7)

            hasdata = .false.
         ELSE
            WRITE(*,*) " Unexpected header ", header
            CALL EXIT(-1)
         ENDIF
      ENDDO
      IF (nat > 0) DEALLOCATE(atoms, forces, msgbuffer)
      END PROGRAM
