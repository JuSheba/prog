module fourier_module
  implicit none
  contains

    recursive subroutine FFT_DIF(x, key_sign)
      integer :: key_sign, N, i
      complex :: a
      complex, dimension(:) :: x
      complex, dimension(:), allocatable :: even, odd
      real(8) :: pi
     !________________________________________________
      pi = 4d0 * atan(1d0)
      N  = size(x)

      if(N <= 1) return

      allocate(even(N/2), odd((N+1)/2))
      even = x(1:N:2)
      odd  = x(2:N:2)

      call FFT_DIF(even, key_sign)
      call FFT_DIF(odd, key_sign)

      do i = 1, N/2
        a = exp(complex(0d0, key_sign*2d0*pi*i) / N) * even(i)
        x(i)     = odd(i) + a
        x(i+N/2) = odd(i) - a
      end do

      deallocate(even)
      deallocate(odd)

    end subroutine FFT_DIF

end module fourier_module
