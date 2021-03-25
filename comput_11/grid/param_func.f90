module param_func
  implicit none
  ! Краевые условия:
  ! t лежит в отрезке [a, b]. Значение x(t) в точках: x(a) = alpha, x(b) = beta
  ! num_int - кол-во интервалов, на которые разбивается [a, b].
  ! Пример задачи - спутник должен пролететь из одной точки в другую за опр.
  ! промежуток времени.


  real(8), parameter :: a = 1d0, b = 2d0
  integer, parameter :: num_int = 100
  real(8), parameter, dimension(1:1) :: alpha = (/0.0d0/), beta = (/0.5d0/)

contains

  function func(t, x, dx)
    implicit none
    real(8), dimension(:) :: x, dx
    real(8), dimension(size(x)) :: func
    real(8) :: t
      func = -2*t*(dx)**2
  end function func

  function check_result(t) result(solution)
    implicit none
    real(8) :: solution, t
      solution = 1d0 - 1d0 / t
    end function check_result

end module param_func
