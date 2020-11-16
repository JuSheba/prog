module IterativeModule
implicit none
contains
!_______________________________________________________________________
! 1. Jacobi method
!_______________________________________________________________________
subroutine Jacobi(n, A, B, C, res,eps)

  integer(8)  :: n, i, score, scoreMax
  real(8), dimension(:), allocatable      :: B, res
  real(8), dimension(1:size(B),1:size(B)) :: D, Z, invertD
  real(8), dimension(:,:), allocatable    :: A, C
  real(8), dimension(1:size(B))           :: G, res0
  real(8) :: eps

  write(*,*) 'You chose Jacobi method.'
  write(*,*) 'Calculation solution vector..'

  do i = 1,n
    D(i,i) = C(i,i)
    invertD(i,i) = 1/D(i,i)
  end do

  Z = matmul(invertD,(D - A))

  G = matmul(invertD, B)

  score = 0
  scoreMax = 100
  res0 = res

  do while (sqrt(sum((res0-res)**2))>eps .and. score<scoreMax)
    res0 = res
    res = matmul(Z,res) + G
    score = score + 1
  end do

end subroutine Jacobi
!_______________________________________________________________________
! 2. Seidel method
!_______________________________________________________________________
subroutine Seidel(n, A, B, C, res, eps)

  integer(8)  :: n, i, score, scoreMax
  real(8), dimension(:,:), allocatable :: A, C
  real(8), dimension(:), allocatable   :: res, B
  real(8), dimension(1:size(B))        :: res0
  real(8) :: eps

  write(*,*) 'You chose Seidel method.'
  write(*,*) 'Calculation solution vector..'

  score    = 0
  scoreMax = 100
  res0 = res + 1d0

  do while (sqrt(sum((res0-res)**2))>eps .and. score<scoreMax)
    res0 = res
    do i = 1,n
        res(i) = (B(i) - dot_product(A(i,1:i-1),res(1:i-1)) &
        - dot_product(A(i,i+1:n),res0(i+1:n))) / A(i,i)
    end do
    score = score + 1
  end do

end subroutine Seidel
!_______________________________________________________________________
! 3. Successive over-relaxation
!_______________________________________________________________________
subroutine SOR(n, A, B, C, res, eps)

  integer(8)  :: n, i, j, score, scoreMax, max
  real(8), dimension(:,:), allocatable :: A, C
  real(8), dimension(:), allocatable   :: res, B
    real(8), dimension(1:size(B),1:size(B)) :: P
  real(8), dimension(1:size(B))        :: res0, Q, Q0, maxima
  real(8) :: eps

  write(*,*) 'You chose SOR method.'
  write(*,*) 'Calculation solution vector..'

  res = 0

  P(i,j) = -A(i,j)/A(i,i)

  Q(i) = B(i)/A(i,i)

  Q0 = Q

  score = 0
  scoreMax = 100


  do i = 1,n
    maxima(i) = abs(Q(i))
  end do
  max = sum(maxloc(maxima))         ! мне очень стыдно за весь этот костыль


  do while (abs(Q(max))>eps .and. score<scoreMax)
    Q0 = Q
    do concurrent(i=1:n)
      Q(i) = Q0(i) + P(i,max)*Q0(max)
    end do
    score = score + 1
  end do





end subroutine SOR
!_______________________________________________________________________
end module IterativeModule
