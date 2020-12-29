program quad_gauss
  use legendre_mod
  use elimination_mod
  !_________________________________________________
  implicit none
    !contains
      integer(8) :: n, i
      real(8) :: eps
      real(8), allocatable :: A(:), t(:)
      real(8), allocatable :: M(:,:)
      character(2) :: num
      !______________________________________________________________
      call getarg(1, num)
      read(num,'(i2)') n

      allocate(A(1:n), t(1:n), M(0:n-1,0:n))
      call root_legendre(t, n)
        
      forall (i = 0:n-1) M(i,0:n-1) = t**i
      M(1::2,n) = 0d0
      forall (i = 0:n-1:2) M(i,n) = 2d0 / (i+1)
      
      eps = epsilon(1d0)
      call GaussJordan(n, eps, M, A)

      open(22,file='quad'//trim(num)//'.dat')
      do i = 1, n
        write(22,*) A(i), t(i)
      end do
      close(22)

end program quad_gauss
