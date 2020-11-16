module newtonModule
  use Elimination
  implicit none
!________________________________________________________________________________________
  contains
      interface
        function func
          integer :: i,n
          real(8), dimension(:) :: X
          real(8), dimension(size(X)) :: func
            n = size(X)
            forall(i = 1:n) func(i) = X(i)**2
        end function func
      end interface

  !______________________________________________________________________________________

  function der(func, x0, n)
    integer :: i, j, n
    real(8) :: x0(n), f_x0(n), f_x1(n), f_x2(n), der(n,n), x1(n), x2(n)

    x1 = x0
    x2 = x0
    t  = 0.00001d0

    do i = 1, n
      x1(i) = x0(i) + t
      x2(i) = x0(i) + 2d0*t
      f_x0 = func(x0)
      f_x1 = func(x1)
      f_x2 = func(x2)
      do j = 1, n
        der(j,i) = -1.5d0 * f_x0(j)/t + &
                        2d0   * f_x1(j)/t - &
                        0.5d0 * f_x2(j)/t
      end do
  end function der
  !______________________________________________________________________________________

  function mNewton(func, x0, max)
    integer :: i, n, max
    real(8) :: x0(:), mNewton(size(x0)), eps
    real(8), allocatable, dimension(:)   :: f_xj_0, xj_0, xj_1
    real(8), allocatable, dimension(:,:) :: C

    eps = epsilon(1d0)
    n = size(x0)
    i = 1

    allocate(C(n,n))
    allocate(f_xj_0(n), xj_0(n), xj_1(n))

    xj_0 = x0
    f_xj_0 = func(xj_0)
    C = der(func, xj_0, n)
    call GaussPivot(n, eps, C, xj_1 )
    xj_1 =

    do while ((maxval(abs(xj_1 - xj_0)) >= eps).and.(i < max))
      xj_0 = xj_1
      f_xj_0 = func(xj_0)
      C = der(func, xj_0, n)
      call GaussPivot(n, eps, C, xj_1 )
      xj_1 =
      i = i+1
    end do
    mNewton = xj_1
  end function mNewton

end module newtonModule