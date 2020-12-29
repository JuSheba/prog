module bernoulli_mod
  use gorner_mod
  implicit none
  contains

    function bernoulli(M)
      real(8)                :: x0, x1, eps
      integer(8)             :: i, n
      real(8), dimension(0:) :: M
      real(8)                :: y(0:size(M)), bernoulli(0:size(M)-2)
      !______________________________________________________________
      eps = epsilon(1d0)
      n   = size(M) - 1
      call random_number(y)

      do i = 0, n - 1
        call random_number(y(0:n-i))
        y(0) = -dot_product(M(1:n-i),y(1:n-i)) / M(0)
        x1 = y(0) / y(1)
        x0 = y(1)
        do while (abs(x1-x0) > eps)
          x0 = x1
          y  = cshift(y(0:n-i),-1)
          y(0) = -dot_product(M(1:n-i),y(1:n-i)) / M(0)
          x1   = y(0)/y(1)
        end do
        bernoulli(i) = x1
        M(0:n-i-1) = gorner(x1, M(0:n-i))
      end do
      !write(*,*)
    end function bernoulli

end module bernoulli_mod
