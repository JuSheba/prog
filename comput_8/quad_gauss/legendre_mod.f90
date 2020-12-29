module legendre_mod
  use bernoulli_mod
  contains

    recursive subroutine coef_legendre(M, n)
    ! Рекурсивная сабрутина для поиска коэффициентов полинома Лежандра заданной степени n
    ! Первые два коэф. задаём явно,(ф-ла работает для n > 1), остальные ищем реккурентно.
    implicit none
    integer(8) :: n, i
    real(8), dimension(0:n)     :: M
    real(8), dimension(0:n-1)   :: poli_1
    real(8), dimension(0:n-2)   :: poli_2

    select case(n)
    case(0)
      M(0) = 1
    case(1)
      M(0) = 1d0
      M(1) = 0d0
    case default
      call coef_legendre(poli_1, n-1)
      call coef_legendre(poli_2, n-2)

      M(0) = (2d0*n-1d0) / n*poli_1(0)
      m(1)=(2d0*n-1d0)  /  n*poli_2(1)

      do i = 2, n -1
         M(i) = (2d0 * n-1d0) / n * poli_1(i) - (n - 1d0) / n * poli_2(i-2)
      end do
      M(n) = -(n - 1d0) / n * poli_2(n-2)

    end select

    end subroutine coef_legendre

!_______________________________________________________________________________________

    subroutine root_legendre(x, n)
    ! Процедура, вычисляющая корни полинома Лежандра. Помним о том, что при нечётной
    ! степени есть нулевой корень, кроме нулевого корня все остальные - парные.
    implicit none
    integer(8) :: n, i
    real(8), dimension(0:n) :: M
    real(8), dimension(0:n-1) :: x
    real(8), dimension(0:n/2-1) :: mod_root

      call coef_legendre(M, n)

      mod_root = sqrt(bernoulli(M(0:2)))

      if (mod(n,2) == 0) then
        x = (/-mod_root, mod_root(n/2-1:0:-1)/)
      else
        x = (/-mod_root, 0d0, mod_root(n/2-1:0:-1)/)
      end if

    end subroutine root_legendre

end module legendre_mod
