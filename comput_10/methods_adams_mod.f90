module methods_adams_mod
  use param_func_mod
  use rk_mod
  use coeff_adams_mod
  use newton_mod
  contains

  function e_adams(x0, t0, x)
    implicit none
    integer :: i, n
    real(8) :: t0
    real(8), dimension(:)   :: x0, e_adams(1:size(x0)), p(1:size(x0))
    real(8), dimension(:,:) :: x(-e_n+1 : 0,1:size(x0))

    n = e_n
    p = 0

    do i = -n+1, 0
      p = p + A(n,-i) * func(t0 + h*i, x(i,:))
    end do

    e_adams = x0 + h*p

  end function e_adams

  function i_adams(x0, t0, x)
    implicit none
    integer :: i, n
    real(8) :: t0
    real(8), dimension(:)   :: x0, i_adams(1:size(x0)), p(1:size(x0))
    real(8), dimension(:,:) :: x(-i_n+2:0, 1:size(x0))

    n = i_n
    p = 0d0

    do i = -n+2, 0
      p = p + B(n,-i) * func(t0 + h*i, x(i,:))
    end do

    call newton(x0, 50, i_adams, fix_func)

  contains

    function fix_func(x)
      implicit none
      real(8), dimension(:) :: x, fix_func(1:size(X))

      fix_func = p*h + h*B(i_n, -1) * func(t0+h, x) + x0 - x
    end function fix_func

  end function i_adams


end module methods_adams_mod
