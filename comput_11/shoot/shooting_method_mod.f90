module shooting_method_mod
  use rk_mod
  use newton_mod
  use param_func
  implicit none

contains

  function shooting(func, a, b, alpha, beta, num_int, key_word)
    implicit none
    interface
      function func(t, x)
        implicit none
        real(8), dimension(:) :: x
        real(8), dimension(size(x)) :: func
        real(8) :: t
      end function func
    end interface
    real(8), dimension(1:) :: alpha, beta, x_res(1:2)
    real(8), dimension(1:size(alpha)) :: x_0, x_k, elem, zero
    integer :: num_int, n
    real(8) :: t_0, t_k, h1, a, b, shooting, t
    character(len=*) :: key_word

    select case(key_word)
    case('left')
      n = size(alpha)
      h1 = (b - a) / num_int
      x_0 = alpha; x_k = beta
      t_0 = a;     t_k = b
    case('right')
      n = size(alpha)
      h1 = - (b - a) / num_int
      x_0 = beta; x_k = alpha
      t_0 = b;    t_k = a
    end select

    zero = 0d0
    elem(1:n) = 0d0
    x_res(1:n) = x_0
    call newton(zero, 250, elem, get_rk)
    x_res(n+1:2*n) = elem
    t = t_0

  open(11, file = 'shoot.dat')
    do while (t <= t_k)
      write(11,*) t, x_res(:)
      x_res = rk(x_res, t, h1)
      t = t + h1
    end do
  close(11)

  contains

    function get_rk(xx_0)
      real(8) :: t
      real(8) :: xx_0(1:)
      real(8) :: get_rk(1:size(xx_0)), x(1:2*size(xx_0))
        x(1:n) = x_0
        x(n+1:2*n) = xx_0
        t = t_0
        do while (t <= t_k)
          x = rk(x, t, h1)
          t = t + h1
        end do
        get_rk = x(1:n) - x_k
    end function get_rk

  end function shooting

end module shooting_method_mod
