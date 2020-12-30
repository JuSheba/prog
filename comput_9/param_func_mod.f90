module param_func_mod
  implicit none
  real(8) :: data_t = 3d0
  real(8), dimension(1:2), parameter :: data = (/ 0d0, 1d0 /)
  real(8) :: h = 2e-4
  integer :: e_n = 5, i_n = 5 ! extra.adams & inter.adams
  contains

  function func(t, x) result(y)
    implicit none
    real(8), dimension(:) :: x, y(1:size(x))
    real(8) :: t

    y(1) =  x(2)
    y(2) =  x(1) + x(2) + cos(t)

end function func

end module param_func_mod
