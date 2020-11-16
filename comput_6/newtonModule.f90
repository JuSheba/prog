module newtonModule
  use Elimination
  use functionModule
  implicit none
!________________________________________________________________________________________
  contains

  function der(func, x0, n)
    integer(8) :: i, j, n
    real(8) :: x0(n), x1(n), der(n,n)
    real(8) :: t
    interface
      function func(X)
        integer :: i, n
        real(8), dimension(:) :: X
        real(8), dimension(size(X)) :: func
      end function func
    end interface

    x1 = x0
    t  = 1e-5

    do i = 1, size(x0)
      x1(i) = x0(i) + t
      der(:,i) = func(x1)/t - func(x0)/t
    end do

  end function der
  !______________________________________________________________________________________

  function mNewton(func, x0, max)
    integer(8) :: i, n, max
    real(8) :: x0(:), x1(size(x0)), mNewton(size(x0))
    real(8) :: eps

    interface
      function func(X)
        integer(8) :: i,n
        real(8), dimension(:) :: X
        real(8), dimension(size(X)) :: func
      end function func
    end interface

    eps = epsilon(1d0)
    n = size(x0)
    i = 0

    x1 = x0

    do while ((maxval(abs(x1 - x0)) >= eps) .and. (i < max))
      x0 = x1
      call Gauss(n, eps, der(func, x0, n), -func(x0), x1)
      x1 = x0 + x1
      i = i + 1
    end do
    mNewton = x1
  end function mNewton

end module newtonModule
