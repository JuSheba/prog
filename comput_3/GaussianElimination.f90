program GaussElimination
use Elimination
implicit none
!_______________________________________________________________________
! DESCRIPTION OF THE PROGRAM :
! This program finds a solution to a linear equations system. The matrix
! of coefficients and the vector of values are read.
! When you start the program, you must select one of the methods for
! solving equations system : Gauss Elimination, Gauss elimination
! with partial pivoting or Gauss-Jordan elimination
! (key : gauss, gausspivot, gaussjordan).
! The result of the program is a vector of solutions.
!_______________________________________________________________________
! Variables used :
!  n - matrix size
!  A - system matrix
!  B - column vector
!  C - augmented system matrix
!  R - residual
!  eps - if elem. in denominator <eps, then program displays a warning
!  res - solution vector
!_______________________________________________________________________
integer(8)  :: n, i, j, k
real(8), dimension(:,:), allocatable :: A, C
real(8), dimension(:), allocatable   :: B, res, R, t
real(8) :: eps, s
character(len=72) :: key
eps = 1e-4
write(*,*) 'You can use : gauss, gausspivot, gaussjordan'
!_______________________________________________________________________
! READING А & В
!_______________________________________________________________________
write(*,*) 'Readind data.dat..'
open(12, file='data.dat')
  read(12, '(2x, i10)')  n
  allocate(A(n,n), B(n))

  do i=1,n
      read(12,*) A(i,1:n)
  enddo

  do i=1,n
    read(12,*) B(i)
  enddo
close(12)

!A = transpose(A)
!_______________________________________________________________________
! GET THE AUGMENTED SYSTEM MATRIX
!_______________________________________________________________________
allocate(C(n,n+1))

do i = 1,n
  do j = 1,n
    C(i,j) =  A(i,j)
  end do
end do

do i=1,n
  C(i,n+1) = B(i)
end do

allocate(res(n))
!_______________________________________________________________________
! CALCULATION SOLUTION
!_______________________________________________________________________
call getarg(1, key)

select case(key)
  case('gauss')
    call Gauss(n,eps,C,res)
  case('gausspivot')
    call GaussPivot(n,eps,C,res)
  case('gaussjordan')
    call GaussJordan(n,eps,C,res)
endselect

open(22, file='result.dat')
write(22,'("#", 1x, i0)') n
do i=1,n
 write(22,*) res(i)
end do

write(*,*) 'Result in result.dat'
allocate(R(n))
R = (sum((matmul(A, res) - B)**2))**0.5
write(*,*) 'Residual R =', R

end program GaussElimination
