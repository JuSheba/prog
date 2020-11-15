program splineFit
  use TDMA
  implicit none

  integer :: n, i, j
  real(8) :: bar, t, h, x_i, f
  real(8), allocatable, dimension(:)   :: X, Y, R, P, S, Q
  real(8), allocatable, dimension(:,:) :: A, B, trB, C, array
  !______________________________________________________________________________________
  ! ЧТЕНИЕ DATA.DAT (X, Y, P)____________________________________________________________
  open(11, file='data.dat')
    read(11, "(2x,i6)") n
    allocate(X(0:n), Y(0:n), R(0:n), P(0:n), S(0:n), Q(0:n))
    allocate(A(0:n,3), B(0:n,3), trB(0:n,3), C(0:n,5), array(0:n,3))
    do i = 0, n
      read(11,*) X(i), Y(i), P(i)
    end do
  close(11)

  ! ЗAПОЛНЕНИЕ МАТРИЦ A, B и Q___________________________________________________________
  do i = 2 , n-2
    A(i,1) = X(i) - X(i-1)           ! верхняя диагональ
    A(i,2) = 2d0 * (X(i+1) - X(i-1)) ! главная диагональ
    A(i,3) = X(i+1) - X(i)           ! нижняя  диагональ
  end do

  A(0,1)   = 2d0 * (X(1) - X(0))
  A(1,2)   = 2d0 * (X(2) - X(0))
  A(1,3)   = X(2) - X(1)
  A(n-1,1) = X(n-1) - X(n-2)
  A(n-1,2) = 2d0 * (X(n) - X(n-2))
  A(n,2) = 2 * (X(n) - X(n-1))

  do i = 1, n-1
    B(i,1) = 1d0  / (X(i)   - X(i-1))                           ! верхняя диагональ
    B(i,2) = -(1d0 / (X(i)   - X(i-1)) + 1d0 / (X(i+1)- X (i))) ! главная диагональ
    B(i,3) = 1d0  / (X(i+1) - X(i))                             ! нижняя диагональ
  end do

  do j = 0, n
    Q(j) = 1d0/P(j)
  end do

  ! ТРАНСПОНИРОВАНИЕ B -> trB_____________________________________________________________
   trB = B
   trB(0,2) = trB(1,1)
   trB(1,1) = 0

   do i = 1, n-2
     bar = trB(i,3)
     trB(i,3) = trB(i+1,1)
     trB(i+1,1) = bar
  enddo

  trB(n,2) = trB(n-1,3)
  trB(n-1,3) = 0

  ! УМНОЖЕНИЕ trB на Q___________________________________________________________________
  do i = 0, n
    do j = 1, 3
      array(i,j) = trB(i,j) / P(i)
    end do
  end do

  ! ЧАСТЬ С УМНОЖЕНИЯМИ МАТРИЦ И ВЕКТОРОВ________________________________________________
  array(n,1) = array(n,2)
  array(n,2) = array(n,3)
  array(n,3) = 0

  C = 6d0 * tridiagmatrixMultiplication(B, array, n)

  do i = 2, n-2
    C(i,2) = C(i,2) + A(i,1)
    C(i,3) = C(i,3) + A(i,2)
    C(i,4) = C(i,4) + A(i,3)
  end do

  C(0,1) = C(0,1) + A(0,1)
  C(0,2) = C(0,2) + A(0,2)
  C(1,1) = C(1,1) + A(1,1)
  C(1,2) = C(1,2) + A(1,2)
  C(1,3) = C(1,3) + A (1,3)
  C(n-1,2) = C(n-1,2) + A(n-1,1)
  C(n-1,3) = C(n-1,3) + A(n-1,2)
  C(n-1,4) = C(n-1,4) + A(n-1,3)
  C(n,3)   = C(n,3)   + A(n,2)

  do i = 1, n-1
    R(i) = B(i,1) * Y(i-1) + B(i,2) * Y(i) + B(i,3) * Y(i+1)
  end do
  R(0) = 0d0
  R(n) = 0d0
  R = 6d0 * R

  C(0,4)   = 0d0
  C(0,5)   = 0d0
  C(1,5)   = 0d0
  C(n-1,5) = 0d0
  C(n,4)   = 0d0
  C(n,5)   = 5d0

  ! ВЕКТОР РЕЗУЛЬТАТОВ R_________________________________________________________________
  S = run(C, R, n)

  do i = 1, n-1
    R(i) = array(i,1) * S(i-1) + array(i,2) * S(i) + array(i,3) * S(i+1)
  end do

  R(0) = array(0,2) * S(1)
  R(n) = array(n,1) * S(n-1)

  R = Y - R

  ! АППРОКСИМАЦИЯ________________________________________________________________________
  open(22, file='result.dat')
  do i = 0, 100*n
    x_i = (X(n) - X(0)) * i &
         /(100d0*n) + X(0)
    j = i / 100d0
    h = X(j+1) - X(j)
    t = (x_i - X(j)) / h
    f = R(j) * (1d0 - t) + R(j+1) * t - &
    h * h * t * (1d0 - t) / &
    6d0*((2-t) * S(j) + (1+t) * S(j+1))

    write(22,*) x_i, f
  enddo
 close(22)

end program splineFit
