subroutine wrap_dpbstf(AB, UPLO, N, LDAB, INFO)
implicit none
double precision, intent(inout) :: AB(LDAB,N)
character, intent(in) :: UPLO
integer, intent(in) :: N, LDAB
integer, intent(out) :: INFO
integer :: KD
KD = LDAB - 1
call dpbstf(UPLO,N,KD,AB,LDAB,INFO)
if (INFO < 0) then
	print *, 'ERROR from DPBSTF: the', -INFO, 'argument had an illegal value'
!else if (INFO > 0) then
!	print *, 'ERROR from DPBSTF: the factorization could not be completed because the update element a(i,i) was negative; the matrix A is not positive definate. i =', INFO
end if
end subroutine wrap_dpbstf

subroutine wrap_dsbgst(AB, BB, UPLO, N, LDAB, LDBB, INFO)
implicit none
double precision, intent(inout) :: AB(LDAB,N)
double precision, intent(inout) :: BB(LDBB,N)
character, intent(in) :: UPLO
integer, intent(in) :: N, LDAB, LDBB
integer, intent(out) :: INFO
double precision :: X, WORK(2*N)
integer :: KA, KB, LDX
KA = LDAB - 1
KB = LDBB - 1
LDX = 1
call dsbgst('N',UPLO,N,KA,KB,AB,LDAB,BB,LDBB,X,LDX,WORK,INFO)
if (INFO < 0) then
	print *, 'ERROR from DSBGST: the', -INFO, 'argument had an illegal value'
end if
end subroutine wrap_dsbgst

subroutine wrap_dgbtrf(M, N, KL, KU, AB, LDAB, IPIV, INFO)
implicit none
integer, intent(in) :: M, N, KL, KU
double precision, intent(inout) :: AB(LDAB,N)
integer, intent(in) :: LDAB
integer, intent(out) :: IPIV(N)
integer, intent(out) :: INFO
call dgbtrf(M, N, KL, KU, AB, LDAB, IPIV, INFO)
if (INFO < 0) then
	print *, 'ERROR from DGBTRF: the', -INFO, 'argument had an illegal value'
end if
end subroutine wrap_dgbtrf

subroutine wrap_dgbtrs(TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO)
implicit none
character*1, intent(in) :: TRANS
integer, intent(in) :: N, KL, KU, NRHS
double precision, intent(in) :: AB(LDAB,N)
integer, intent(in) :: IPIV(N)
double precision, intent(inout) :: B(LDB,NRHS)
integer, intent(in) :: LDAB, LDB
integer, intent(out) :: INFO
call dgbtrs(TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO)
if (INFO < 0) then
	print *, 'ERROR from DGBTRS: the', -INFO, 'argument had an illegal value'
end if
end subroutine wrap_dgbtrs
