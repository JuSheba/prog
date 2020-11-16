program MethodIterative
use IterativeModule
implicit none
!_______________________________________________________________________
! DESCRIPTION OF THE PROGRAM :
!_______________________________________________________________________
! Variables used :
!  n - matrix size
!  A - system matrix
!  B - column vector
!  C - augmented system matrix
!  R - residual
!_______________________________________________________________________
integer(8) :: n, i, j, bar, bar2
real(8)    :: eps
real(8), dimension(:,:), allocatable :: A, C
real(8), dimension(:), allocatable   :: B, res, R

character(len=72) :: key

eps = 1e-4
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
! TEST DIAGONALLY DOMINANT MATRIX
!_______________________________________________________________________
write(*,*) 'Testing diagonally dominant matrix..'
bar2 = 0
do i = 1,n
  bar = sum(abs(C(i,:)))
  bar = bar - abs(C(i,i))
  if(abs(C(i,i)) <= bar) then
    bar2 = 1
  end if
  bar = 0
end do

if (bar == 1) then
  write(*,*) 'This matrix does not have a diagonally dominant. &
   Program exited.'
  call exit()
end if



write(*,*) 'This matrix has a diagonally dominant.'

!_______________________________________________________________________
! CALCULATION SOLUTION
!_______________________________________________________________________
call getarg(1, key)

select case(key)
case('jacobi')
!    allocate(D(n,n), invertD(n,n), Z(n,n))
!    allocate(G(n))
    call Jacobi(n, A, B, C, res,eps)
  case('seidel')
    call Seidel(n, A, B, C, res, eps)
  case('sor')
    call SOR(n, A, B, C, res, eps)
endselect
!_______________________________________________________________________
! DATA OUTPUT
!_______________________________________________________________________
open(22, file='result.dat')
write(22,'("#", 1x, i0)') n
do i=1,n
 write(22,*) res(i)
end do

write(*,*) 'Result in result.dat'
allocate(R(n))
R = (sum((matmul(A, res) - B)**2))**0.5
write(*,*) 'Residual R =', R

end program MethodIterative
