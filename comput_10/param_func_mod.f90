module param_func_mod
  implicit none
  real(8), parameter :: t_0 = 0d0, t_k = 1d2, h = 3e-3
  integer, parameter :: e_n = 3, i_n = 3
  real(8), dimension(1:2), parameter :: X0 = (/ -3d-1, 0d-1 /)

contains
    function func(t, x)
      implicit none
        real(8) :: t
        real(8), parameter :: k = 0.28d0
        real(8), dimension(:) :: x
        real(8), dimension(1:size(x)) :: func
          func(1) = x(2)
          func(2) = k * (1d0 - x(1)**2) * x(2) - x(1)
    end function func

end module param_func_mod
