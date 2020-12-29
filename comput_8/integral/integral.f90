program integral
  use func_mod
  implicit none
  real(8) :: a = 1.0, b = 3.0
  real(8) :: result
  integer :: n = 6
    call gauss_integral(a, b, n, func, result)
end program integral

    subroutine gauss_integral(a0, b0, n, func, res)
      implicit none
      real(8) :: a0, b0, res
      integer :: n, i
      real(8), dimension(1:n) :: A, t
      character(2) :: num
      interface
        function func(x)
        real(8) :: x
        real(8) :: func
        end function
      end interface

        open(11,file='quad6.dat')
          do i = 1, n
            read(11,*) A(i), t(i)
          end do
        close(11)

        res = 0
        do i = 1, n
            res = res + A(i) * (b0 - a0) / 2 * func(t(i)*(b0-a0)/2 + (a0+b0)/2)
        end do

        write(*,*) res
    end subroutine gauss_integral
