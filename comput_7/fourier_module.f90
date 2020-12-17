module fourier_module
  implicit none
  contains

    recursive subroutine FFT_DIF(x, y, w)
      integer :: N, i
      complex, dimension(0:) :: x, w
      complex, dimension(0:size(x)-1) :: y
     !________________________________________________
      N  = size(x)

      if(N > 2) then
        call FFT_DIF(x(0:N/2-1) + x(N/2:N-1),&
                     y(0:N-2:2), w(0:N-1:2))
        call FFT_DIF(w(0:N/2-1)*(x(0:N/2-1) - x(N/2:N-1)),&
                     y(1:N-1:2), w(0:N-1:2))
      else
        do i = 0, 1
          y(i) = x(0) + (-1)**i*x(1)
        end do
      end if

    end subroutine FFT_DIF

end module fourier_module
