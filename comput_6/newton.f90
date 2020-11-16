program newton
  use functionModule
  use newtonModule
  implicit none

  real(8),allocatable,dimension(:) :: x0, x
  integer(8) :: i, max, n

  n = 200 ! size x
  max = 10e6

  allocate(x(n), x0(n))
  forall(i = 1:n) x0(i) = 0.5

  x = mNewton(func, x0, max)
  write(*,*)  sqrt(dot_product(func(x), func(x)))
  open(12, file='result.dat')
    do i = 1, n
      write(12,*) x
    end do
  close(12)

end program newton
