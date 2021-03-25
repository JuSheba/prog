module param_func
  implicit none
  ! Краевые условия:
  ! t лежит в отрезке [a, b]. Значение x(t) в точках: x(a) = alpha, x(b) = beta
  ! num_int - кол-во интервалов, на которые разбивается [a, b].
  ! Пример задачи - спутник должен пролететь из одной точки в другую за опр.
  ! промежуток времени.


  real(8), parameter :: a = 1d0, b = 2d0
  integer, parameter :: num_int = 20
  real(8), parameter, dimension(1:1) :: alpha = (/0.7d0/), beta = (/0.5d0/)

contains

  function func(t, x)
    implicit none
    real(8), dimension(:) :: x
    real(8), dimension(size(x)) :: func
    real(8) :: t
      func(1) = x(2)
      func(2) = -2*t*(x(2))**2
  end function func

end module param_func
