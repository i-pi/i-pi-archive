! The main program which runs our driver test case potentials
! 
! Copyright (C) 2013, Joshua More and Michele Ceriotti
! 
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
! 
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
!
! Currently the potentials implemented are the Lennard-Jones
! potential, the Silvera-Goldman para-hydrogen potential and
! the ideal gas (i.e. no interaction at all)

      PROGRAM DRIVER
         USE LJ
         USE SG
         USE F90SOCKETS, ONLY : open_socket, writebuffer, readbuffer
      IMPLICIT NONE

      ! SOCKET VARIABLES
      INTEGER, PARAMETER :: MSGLEN=12   ! length of the headers of the driver/wrapper communication protocol
      INTEGER socket, inet, port        ! socket ID & address of the server
      CHARACTER(LEN=1024) :: host

      ! COMMAND LINE PARSING
      CHARACTER(LEN=1024) :: cmdbuffer
      INTEGER ccmd, vstyle
      LOGICAL verbose
      INTEGER commas(4), par_count      ! stores the index of commas in the parameter string
      DOUBLE PRECISION vpars(4)         ! array to store the parameters of the potential

      ! SOCKET COMMUNICATION BUFFERS
      CHARACTER(LEN=12) :: header
      LOGICAL :: isinit=.false., hasdata=.false.
      INTEGER cbuf
      CHARACTER(LEN=2048) :: initbuffer      ! it's unlikely a string this large will ever be passed...
      DOUBLE PRECISION, ALLOCATABLE :: msgbuffer(:)

      ! PARAMETERS OF THE SYSTEM (CELL, ATOM POSITIONS, ...)
      DOUBLE PRECISION sigma, eps, rc, rn, ks ! potential parameters
      INTEGER nat
      DOUBLE PRECISION pot, dpot
      DOUBLE PRECISION, ALLOCATABLE :: atoms(:,:), forces(:,:), datoms(:,:)
      DOUBLE PRECISION cell_h(3,3), cell_ih(3,3), virial(3,3), mtxbuf(9), dip(3)
      DOUBLE PRECISION volume
      DOUBLE PRECISION, PARAMETER :: fddx = 1.0d-5

      ! NEIGHBOUR LIST ARRAYS
      INTEGER, DIMENSION(:), ALLOCATABLE :: n_list, index_list
      DOUBLE PRECISION init_volume, init_rc ! needed to correctly adjust the cut-off radius for variable cell dynamics
      DOUBLE PRECISION, ALLOCATABLE :: last_atoms(:,:) ! Holds the positions when the neighbour list is created
      DOUBLE PRECISION displacement ! Tracks how far each atom has moved since the last call of nearest_neighbours

      INTEGER i, j 

      ! parse the command line parameters
      ! intialize defaults
      ccmd = 0
      inet = 1
      host = "localhost"//achar(0)
      port = 31415
      verbose = .false.
      par_count = 0
      vstyle = -1 
      rc = 0.0d0
      init_rc = 0.0d0
      volume = 0.0d0
      init_volume = 0.0d0

      DO i = 1, IARGC()
         CALL GETARG(i, cmdbuffer)
         IF (cmdbuffer == "-u") THEN ! flag for unix socket
            inet = 0
            ccmd = 0
         ELSEIF (cmdbuffer == "-h") THEN ! read the hostname
            ccmd = 1
         ELSEIF (cmdbuffer == "-p") THEN ! reads the port number
            ccmd = 2
         ELSEIF (cmdbuffer == "-m") THEN ! reads the style of the potential function
            ccmd = 3
         ELSEIF (cmdbuffer == "-o") THEN ! reads the parameters
            ccmd = 4
         ELSEIF (cmdbuffer == "-v") THEN ! flag for verbose standard output
            verbose = .true.
         ELSE
            IF (ccmd == 0) THEN
               WRITE(*,*) " Unrecognized command line argument", ccmd
               CALL helpmessage
               STOP "ENDED"
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
               ELSEIF (trim(cmdbuffer) == "harm") THEN
                  vstyle = 3
               ELSEIF (trim(cmdbuffer) == "morse") THEN
                  vstyle = 4
               ELSEIF (trim(cmdbuffer) == "zundel") THEN
                  vstyle = 5
               ELSEIF (trim(cmdbuffer) == "qtip4pf") THEN
                  vstyle = 6
               ELSEIF (trim(cmdbuffer) == "gas") THEN
                  vstyle = 0  ! ideal gas
               ELSE
                  WRITE(*,*) " Unrecognized potential type ", trim(cmdbuffer)
                  WRITE(*,*) " Use -m [gas|lj|sg|harm|morse|zundel|qtip4pf] "
                  STOP "ENDED"
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
            ccmd = 0
         ENDIF
      ENDDO

      IF (vstyle == -1) THEN
         WRITE(*,*) " Error, type of potential not specified."
         CALL helpmessage
         STOP "ENDED"
      ELSEIF (0 == vstyle ) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error: no initialization string needed for ideal gas."
            STOP "ENDED" 
         ENDIF   
         isinit = .true.
      ELSEIF (6 == vstyle ) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error:  no initialization string needed for qtip4pf."
            STOP "ENDED" 
         ENDIF 
         isinit = .true.
      ELSEIF (5 == vstyle ) THEN
         IF (par_count /= 0) THEN
            WRITE(*,*) "Error: no initialization string needed for zundel."
            STOP "ENDED" 
         ENDIF   
         CALL prezundelpot()
         CALL prezundeldip()
         isinit = .true.
      ELSEIF (4 == vstyle ) THEN
         IF (par_count == 0) THEN ! defaults (OH stretch)
            vpars(1) = 1.8323926 ! r0
            vpars(2) = 0.18748511263179304 ! D
            vpars(3) = 1.1562696428501682 ! a
         ELSEIF ( 2/= par_count) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For morse potential use -o r0,D,a (in a.u.) "
            STOP "ENDED" 
         ENDIF 
         isinit = .true.
      ELSEIF (vstyle == 1) THEN
         IF (par_count /= 3) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For LJ potential use -o sigma,epsilon,cutoff "
            STOP "ENDED" ! Note that if initialization from the wrapper is implemented this exit should be removed.
         ENDIF   
         sigma = vpars(1)
         eps = vpars(2)
         rc = vpars(3)
         rn = rc*1.2
         isinit = .true.
      ELSEIF (vstyle == 2) THEN
         IF (par_count /= 1) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For SG potential use -o cutoff "
            STOP "ENDED" ! Note that if initialization from the wrapper is implemented this exit should be removed.
         ENDIF
         rc = vpars(1)
         rn = rc*1.2
         isinit = .true.
      ELSEIF (vstyle == 3) THEN
         IF (par_count /= 1) THEN
            WRITE(*,*) "Error: parameters not initialized correctly."
            WRITE(*,*) "For 1D harmonic potential use -o k "
            STOP "ENDED" ! Note that if initialization from the wrapper is implemented this exit should be removed.
         ENDIF
         ks = vpars(1)
         isinit = .true.
      ENDIF

      IF (verbose) THEN
         WRITE(*,*) " DRIVER - Connecting to host ", trim(host)
         IF (inet > 0) THEN
            WRITE(*,*) " on port ", port, " using an internet socket."
         ELSE
            WRITE(*,*) " using an UNIX socket."
         ENDIF
      ENDIF

      ! Calls the interface to the POSIX sockets library to open a communication channel
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
            CALL readbuffer(socket, cbuf)
            CALL readbuffer(socket, initbuffer, cbuf)
            IF (verbose) WRITE(*,*) " Initializing system from wrapper, using ", trim(initbuffer)
            isinit=.true. ! We actually do nothing with this string, thanks anyway. Could be used to pass some information (e.g. the input parameters, or the index of the replica, from the driver
         ELSEIF (trim(header) == "POSDATA") THEN  ! The driver is sending the positions of the atoms. Here is where we do the calculation!

            ! Parses the flow of data from the socket
            CALL readbuffer(socket, mtxbuf, 9)  ! Cell matrix
            cell_h = RESHAPE(mtxbuf, (/3,3/))
            CALL readbuffer(socket, mtxbuf, 9)  ! Inverse of the cell matrix (so we don't have to invert it every time here)
            cell_ih = RESHAPE(mtxbuf, (/3,3/))
            
            ! The wrapper uses atomic units for everything, and row major storage.
            ! At this stage one should take care that everything is converted in the
            ! units and storage mode used in the driver.
            cell_h = transpose(cell_h)
            cell_ih = transpose(cell_ih)
            ! We assume an upper triangular cell-vector matrix
            volume = cell_h(1,1)*cell_h(2,2)*cell_h(3,3)

            CALL readbuffer(socket, cbuf)       ! The number of atoms in the cell
            IF (nat < 0) THEN  ! Assumes that the number of atoms does not change throughout a simulation, so only does this once
               nat = cbuf
               IF (verbose) WRITE(*,*) " Allocating buffer and data arrays, with ", nat, " atoms"
               ALLOCATE(msgbuffer(3*nat))
               ALLOCATE(atoms(nat,3), datoms(nat,3))
               ALLOCATE(forces(nat,3))
               atoms = 0.0d0
               datoms = 0.0d0
               forces = 0.0d0
               msgbuffer = 0.0d0
            ENDIF

            CALL readbuffer(socket, msgbuffer, nat*3)
            DO i = 1, nat
               atoms(i,:) = msgbuffer(3*(i-1)+1:3*i)
            ENDDO

            IF (vstyle == 0) THEN   ! ideal gas, so no calculation done
               pot = 0
               forces = 0.0d0
               virial = 0.0d0
            ELSEIF (vstyle == 3) THEN ! 1D harmonic potential, so only uses the first position variable
               pot = 0.5*ks*atoms(1,1)**2
               forces = 0
               forces(1,1) = -ks*atoms(1,1)
               virial = 0
               virial(1,1) = forces(1,1)*atoms(1,1)
            ELSEIF (vstyle == 4) THEN ! Morse potential. 
               IF (nat/=1) THEN
                  WRITE(*,*) "Expecting 1 atom for 3D Morse (use the effective mass for the atom mass to get proper frequency!) "
                  STOP "ENDED"
               ENDIF
               CALL getmorse(vpars(1), vpars(2), vpars(3), atoms, pot, forces)               
            ELSEIF (vstyle == 5) THEN ! Zundel potential. 
               IF (nat/=7) THEN
                  WRITE(*,*) "Expecting 7 atoms for Zundel potential, O O H H H H H "
                  STOP "ENDED"
               ENDIF
               
               CALL zundelpot(pot,atoms)
               CALL zundeldip(dip,atoms)

               datoms=atoms
               DO i=1,7  ! forces by finite differences
                  DO j=1,3                     
                     datoms(i,j)=atoms(i,j)+fddx
                     CALL zundelpot(dpot, datoms)
                     datoms(i,j)=atoms(i,j)-fddx                     
                     CALL zundelpot(forces(i,j), datoms)
                     datoms(i,j)=atoms(i,j)
                     forces(i,j)=(forces(i,j)-dpot)/(2*fddx)
                  ENDDO
               ENDDO
               ! do not compute the virial term
            ELSEIF (vstyle == 6) THEN ! qtip4pf potential.             
               IF (mod(nat,3)/=0) THEN
                  WRITE(*,*) " Expecting water molecules O H H O H H O H H but got ", nat, "atoms"
                  STOP "ENDED"
               ENDIF
               vpars(1) = cell_h(1,1)
               vpars(2) = cell_h(2,2)
               vpars(3) = cell_h(3,3)
               IF (cell_h(1,2).gt.1d-10 .or. cell_h(1,3).gt.1d-12  .or. cell_h(2,3).gt.1d-12) THEN
                  WRITE(*,*) " qtip4pf PES only works with orthorhombic cells"
                  STOP "ENDED" 
               ENDIF
               CALL qtip4pf(vpars(1:3),atoms,nat,forces,pot,virial)
               
               ! do not compute the virial term
            ELSE
               IF ((allocated(n_list) .neqv. .true.)) THEN
                  IF (verbose) WRITE(*,*) " Allocating neighbour lists."
                  ALLOCATE(n_list(nat*(nat-1)/2))
                  ALLOCATE(index_list(nat))
                  ALLOCATE(last_atoms(nat,3))
                  last_atoms = 0.0d0
                  CALL nearest_neighbours(rn, nat, atoms, cell_h, cell_ih, index_list, n_list)
                  last_atoms = atoms
                  init_volume = volume
                  init_rc = rc
               ENDIF

               ! Checking to see if we need to re-calculate the neighbour list
               rc = init_rc*(volume/init_volume)**(1.0/3.0)
               DO i = 1, nat
                  CALL separation(cell_h, cell_ih, atoms(i,:), last_atoms(i,:), displacement)
                  ! Note that displacement is the square of the distance moved by atom i since the last time the neighbour list was created.
                  IF (4*displacement > (rn-rc)*(rn-rc)) THEN
                     IF (verbose) WRITE(*,*) " Recalculating neighbour lists"
                     CALL nearest_neighbours(rn, nat, atoms, cell_h, cell_ih, index_list, n_list)
                     last_atoms = atoms
                     rn = 1.2*rc
                     EXIT
                  ENDIF
               ENDDO

               IF (vstyle == 1) THEN
                  CALL LJ_getall(rc, sigma, eps, nat, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)
               ELSEIF (vstyle == 2) THEN
                  CALL SG_getall(rc, nat, atoms, cell_h, cell_ih, index_list, n_list, pot, forces, virial)
               ENDIF
               IF (verbose) WRITE(*,*) " Calculated energy is ", pot
            ENDIF
            hasdata = .true. ! Signal that we have data ready to be passed back to the wrapper
         ELSEIF (trim(header) == "GETFORCE") THEN  ! The driver calculation is finished, it's time to send the results back to the wrapper

            ! Data must be re-formatted (and units converted) in the units and shapes used in the wrapper
            DO i = 1, nat
               msgbuffer(3*(i-1)+1:3*i) = forces(i,:)
            ENDDO
            virial = transpose(virial)

            CALL writebuffer(socket,"FORCEREADY  ",MSGLEN)
            CALL writebuffer(socket,pot)  ! Writing the potential
            CALL writebuffer(socket,nat)  ! Writing the number of atoms
            CALL writebuffer(socket,msgbuffer,3*nat) ! Writing the forces
            CALL writebuffer(socket,reshape(virial,(/9/)),9)  ! Writing the virial tensor, NOT divided by the volume
            IF (vstyle==5) THEN ! returns the dipole 
               initbuffer = " "
               WRITE(initbuffer,*) dip(1:3)
               cbuf = LEN_TRIM(initbuffer)
               CALL writebuffer(socket,cbuf) ! Writes back the molecular dipole 
               CALL writebuffer(socket,initbuffer,cbuf)
            ELSE
               cbuf = 7 ! Size of the "extras" string
               CALL writebuffer(socket,cbuf) ! This would write out the "extras" string, but in this case we only use a dummy string.
               CALL writebuffer(socket,"nothing",7)
            ENDIF
            hasdata = .false.
         ELSE
            WRITE(*,*) " Unexpected header ", header
            STOP "ENDED"
         ENDIF
      ENDDO
      IF (nat > 0) DEALLOCATE(atoms, forces, msgbuffer)
 
      CONTAINS
      SUBROUTINE helpmessage
         ! Help banner
         WRITE(*,*) " SYNTAX: driver.x [-u] -h hostname -p port -m [gas|lj|sg|harm|morse|zundel|qtip4pf] "
         WRITE(*,*) "         -o 'comma_separated_parameters' [-v] "
         WRITE(*,*) ""
         WRITE(*,*) " For LJ potential use -o sigma,epsilon,cutoff "
         WRITE(*,*) " For SG potential use -o cutoff "
         WRITE(*,*) " For 1D harmonic oscillator use -o k "
         WRITE(*,*) " For 1D morse oscillator use -o r0,D,a"         
         WRITE(*,*) " For the ideal gas, qtip4pf or zundel no options needed! "
      END SUBROUTINE
   END PROGRAM

    
