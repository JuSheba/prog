module rk_mod
use param_shoot
contains

  function rk(x0, t, h)
    implicit none
    real(8), dimension(:) :: x0
    real(8) :: t, h
    real(8), dimension(1:size(x0)) :: rk, k1, k2, k3, k4


    k1 = h*func(t, x0)
    k2 = h*func(t+h/2d0, x0+K1/2d0)
    k3 = h*func(t+h/2d0, x0+K2/2d0)
    k4 = h*func(t+h,     x0+K3)

    rk = x0 + 1d0 / 6d0 * (k1 + 2d0*k2 + 2d0*k3 + k4)
  end function rk

end module rk_mod
