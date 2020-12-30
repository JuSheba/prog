module newton_mod
use elimination_mod
contains

subroutine newton(x0, count_max, x, func)
  implicit none
  real(8), dimension(1:) :: x0
  real(8), dimension(1:size(x0)) :: x, res
  real(8) :: eps, rand
  integer :: count_max, counter, n
  interface
    function func(x)
      real(8), dimension(1:) :: x
      real(8), dimension(1:size(x)) :: func
  end function func
  end interface
    eps = epsilon(1d0)
    n = size(x0)
    x = x0
    counter = 0
    call random_number(rand)
    res = x0 + rand

    do while (sum(abs(x-res)) > eps .and. counter <= count_max)
        x = res
        call solve(x, res, func)
        counter = counter + 1
    end do
    x = res
end subroutine newton

subroutine solve(x, res, func)
  implicit none
  real(8), dimension(:) :: x
  real(8), dimension(1:size(x)) :: buf
  real(8), dimension(1:size(x)) :: res
  real(8), dimension(1:size(x),1:size(x)) :: jacobi
  real(8), dimension(1:size(x)+1,1:size(x)) :: M
  real(8) :: eps
  integer(8) :: i, n
  interface
    function func(x)
      real(8), dimension(1:) :: x
      real(8), dimension(1:size(x)) :: func
    end function func
  end interface
    n = size(x)
    eps = epsilon(1d0)
    !allocate(buf)

    do i = 1, n
        buf = x
        buf(i) = x(i) + eps
        jacobi(i,1:n) = (func(buf) - func(x)) / eps
    end do

    M(n+1,1:n) = -func(x)
    M(1:n,1:n) = jacobi(1:n,1:n)
    call Gauss(n, eps, M, res)
    res = res + x
end subroutine solve

end module newton_mod
