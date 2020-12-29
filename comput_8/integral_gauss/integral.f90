module gauss_integral_mod
  contains

  subroutine gauss_integral(a,b,n,func,res)

  implicit none
  real(8) :: a, b
  integer :: n, i
  real(8) :: res
  real(8), dimension(1:n) :: A, t
  character(2) :: num
  interface
  function func(x)
  real(8) :: x
  real(8) :: f
  end function
  end interface

  !---------------------------------------
  open(11,file='quad'//trim(num)//'.dat')
    forall(i = 1:n) read(11,*) A(i), t(i)
  close(11)
  !---------------------------------------
  res = 0
  do i = 1, n
      res = res + A(i) * (b - a) / 2 * func(t(i)*(b-a)/2 + (a+b)/2)
  end do

  end subroutine gauss_integral
