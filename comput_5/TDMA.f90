module TDMA
  implicit none
  contains

  ! ФУНКЦИЯ УМНОЖЕНИЯ МАТРИЦ_____________________________________________________________
  function tridiagmatrixMultiplication(A, B, n) result(C)
    integer(4) :: n, i
    real(8)    :: A(0:n,3), B(0:n,3), C(0:n,5)

    do i = 2, n-2
      ! Главная диагональ
      C(i,3) = A(i,1) * B(i-1,3) + A(i,2)*B(i,2) + A(i,3)*B(i+1,1)
      ! Вторая верхняя диагональ
      C(i,2) = A(i,1) * B(i-1,2) + A(i,2)*B(i,1)
      ! Вторая нижняя диагональ
      C(i,4) = A(i,2) * B(i,3) + A(i,3) * B(i+1,2)
      ! Третья верхняя диагональ
      C(i,1) = A(i,1) * B(i-1,1)
      ! Третья нижняя диагональ
      C(i,5) = A(i,3) * B(i+1,3)
    enddo

    C(0,1)   =   A(0,1) * B(0,1) + A(0,2) * B(1,1)
    C(0,2)   =   A(0,1) * B(0,2) + A(0,2) * B(1,2)
    C(0,3)   =   A(0,2) * B(1,3)
    C(1,1)   =   A(1,1) * B(0,1) + A(1,2) * B(1,1)
    C(1,2)   =   A(1,1) * B(0,2) + A(1,2) * B(1,2) + A(1,3) * B(2,1)
    C(1,3)   =   A(1,2) * B(1,3) + A(1,3) * B(2,2)
    C(1,4)   =   A(1,3) * B(2,3)
    C(n-1,1) = A(n-1,1) * B(n-2,1)
    C(n-1,2) = A(n-1,1) * B(n-2,2) + A(n-1,2) * B(n-1,1)
    C(n-1,3) = A(n-1,1) * B(n-2,3) + A(n-1,2) * B(n-1,2) + A(n-1,3) * B(n,1)
    C(n-1,4) = A(n-1,2) * B(n-1,3) + A(n-1,3) * B(n,2)
    C(n,1)   =   A(n,1) * B(n-1,1)
    C(n,2)   =   A(n,1) * B(n-1,2) + A(n,2) * B(n,1)
    C(n,3)   =   A(n,1) * B(n-1,3) + A(n,2) * B(n,2)

  end function tridiagmatrixMultiplication

  ! ФУНКЦИЯ ПРОГОНКИ_____________________________________________________________________
  function run(C, R ,n)
    implicit none
    integer(4) :: n, i
    real(8)    ::  C(0:n,5), R(0:n), alpha(0:n), beta(0:n),&
                   pi(0:n), qi(0:n), ri(0:n), run(0:n)

    beta(0) = 0
    alpha(0) = C(0,1)
    pi(0) = C(0,2)/alpha(0)
    qi(0) = C(0,3)/alpha(0)
    ri(0) = R(0)/alpha(0)

    beta(1) = C(0,2)
    alpha(1)=C(1,2)-pi(0)*beta(1)
    pi(1)=(C(1,3)-qi(0)*beta(1))/alpha(1)
    qi(1)=C(1,4)/alpha(1)
    ri(1)=(R(1)-ri(0)*beta(1))/alpha(1)

    do i = 2, n
      beta(i) = C(i,2)-pi(i-2)*C(i,1)
      alpha(i) = C(i,3)-pi(i-1)*beta(i)-qi(i-2)*C(i,1)
      pi(i) = (C(i,4)-qi(i-1)*beta(i))/alpha(i)
      qi(i) = C(i,5)/alpha(i)
      ri(i) = (R(i)-ri(i-1)*beta(i)-ri(i-2)*C(i,1))/alpha(i)
    end do

    run(n) = ri(n)
    run(n-1) = ri(n-1) - pi(n-1) * run(n)
    do i = n-2, 0, -1
      run(i) = ri(i) - pi(i) * run(i+1) - qi(i) * run(i+2)
    end do
end function run
!________________________________________________________________________________________
end module TDMA
