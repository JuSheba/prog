program newton
  use newtonModule
  implicit none

  real(8),allocatable,dimension(:) :: x0, x
  integer :: i, max, n

  n = 150 ! size x
  max = 10e4

  allocate(x(n), x0(n))
  forall(i = 1:n) x0(i) = 0.5d0

  x = mNewton(func, x0, max)
  write(*,*)  sqrt(dot_product(func(x), func(x)))
  open(12, file='result.dat')
    do i = 1, n
      write(12,*) x
    end do
  close(12)

end program newton
