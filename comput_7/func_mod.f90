module func_mod
  implicit none
  contains

  !real(8) :: pi
!  pi = 4d0 * atan(1d0)

  pure function w_minus(q,p) result(w)
    implicit none
    integer, intent(in) :: q, p
    complex :: w
      w = exp(-2d0 * (4d0 * atan(1d0)) * cmplx(0,1) * q/p)
  end function w_minus

  pure function w_plus(q,p) result(w)
    implicit none
    integer, intent(in) :: q, p
    complex :: w
      w = exp(2d0 * (4d0 * atan(1d0)) * cmplx(0,1) * q/p)
  end function w_plus

end module func_mod
