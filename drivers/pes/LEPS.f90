! Atoms must be in the order A,B,C

SUBROUTINE LEPS_M1(na,atoms,V,F)
  IMPLICIT NONE
  REAL(8),PARAMETER :: abc(3) = (/ .05D0, .05D0, .3D0 /), &
       & d(3) = (/4.746D0, 4.746D0, 3.445D0 /), r_0=.742D0, alpha=1.942D0
  REAL(8), INTENT(IN) :: atoms(3,3) ! A=> atoms(1), B=> atoms(2) C=> atoms(3)
  REAL(8), INTENT(OUT) :: V,F(3,3)
  REAL(8) :: r, Jiii(3)=0D0, J_temp=0D0
  INTEGER :: na, i, ii

  V = 0.0D0
  J_temp = 0.0D0
  do i = 1,3
     ii = mod(i,3)+1
     V = V + Q(atoms(i,:),atoms(ii,:), d(i)) / (1 + abc(i))
     Jiii(i) = J(atoms(i,:),atoms(ii,:), d(i)) / (1 + abc(i))
  end do
  do i = 1,3
     ii = mod(i,3)+1
     J_temp = J_temp + Jiii(i)**2 - Jiii(i)*Jiii(ii)
  end do
  
  V = V - J_temp**0.5
  F =0.0
  
  
CONTAINS
  REAL FUNCTION Q(atom1, atom2, d)
    REAL(8), INTENT(IN) :: d, atom1(3), atom2(3)
    REAL(8) :: r, dist
    r = dist(atom1,atom2)
    Q = d/2 * (3/2 * exp(-2 * alpha * (r-r_0)) - exp(-1 * alpha * (r-r_0)))
  END FUNCTION Q

  REAL FUNCTION J(atom1, atom2, d)
    REAL(8), INTENT(IN) :: atom1(3), atom2(3), d
    REAL(8):: r, dist
    r = dist(atom1,atom2)
    J = d/4D0 * (exp(-2D0 * alpha * (r-r_0)) - 6D0 * exp(-1D0 * alpha * (r-r_0)))
!    write(*,*) atom1(1), atom2(1), d, J

  END FUNCTION J
  
  
  
END SUBROUTINE LEPS_M1


PROGRAM test
  IMPLICIT NONE
  real(8), DIMENSION(3,3) :: atoms, step1, step2, atoms_orig
  real(8) :: V, F(3,3), dist
  integer :: na = 3, i, j
  
  atoms = reshape((/ 1.D0, 1.D0, 1.D0, .0D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /),(/ 3,3 /)) * .0D0
  atoms_orig = atoms
  step1 = reshape((/ 0.D0, 1D0, 1D0, .0D0, 0.D0, 0.D0, .0D0, 0.D0, 0.D0 /),(/ 3,3 /)) * 4D-2
  step2 = reshape((/ 0.D0, 0D0, 1D0, .0D0, 0.D0, 0.D0, .0D0, 0.D0, 0.D0 /),(/ 3,3 /)) * 4D-2
  
  ! atoms = reshape((/ 1.D0, 2.D0, 3.D0, .0D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /),(/ 3,3 /))
  ! CALL LEPS_M1_V(na,atoms,V,F)

  do i=1,101
     do j = 1,101
        CALL LEPS_M1(na,atoms,V,F)
        ! V = 0D0
        ! write(*,*) atoms(1,1), atoms(2,1), atoms(3,1)
        ! write(*,*) dist(atoms(3,:), atoms(2,:)), dist(atoms(1,:), atoms(2,:)), dist(atoms(1,:),atoms(3,:))
        write(*,*) dist(atoms(3,:), atoms(2,:)), dist(atoms(1,:), atoms(2,:)), V
        atoms = atoms + step1
     end do
     atoms = atoms_orig + DBLE(i)*step2
     write(*,*)
  end do

1000 format('3F12.6')
  
END PROGRAM test

REAL(8) FUNCTION dist(p1,p2)
  IMPLICIT NONE
  REAL(8),INTENT(IN),DIMENSION(3) :: p1, p2

  dist = ((p1(1) - p2(1))**2 + (p1(2) - p2(2))**2 + (p1(3) - p2(3))**2D0)**0.5D0

  
END FUNCTION dist
