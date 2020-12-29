module func_mod
  contains
  function func(x)
    implicit none
    real(8), intent(in) :: x
    real(8) :: func
    integer(8) :: i

    func = x**2

  end function func

end module func_mod
