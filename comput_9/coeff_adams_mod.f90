module coeff_adams_mod
  contains

  function A(n, j)   !Часть функции А_nj без рассчёта интеграла, подыинтегральная ф-я
    implicit none      !называется intA, процедура для вычисления интеграла - gauss_integral
    integer :: n, j, i !из 8 домашки и выписанные коэф., чтобы не усложнять себе жизнь.
    real(8) :: A, integral ! Аналогично для B. Спасибо Фортрану за существование contains
    call gauss_integral(0d0, 1d0, intA, integral)

    A = (-1d0)**j / &
        product((/(i, i=1,j)/)) / product((/(i, i=1,n-1-j)/)) * integral
    contains
    function intA(z)
      implicit none
      real(8) :: z, intA
      integer :: i
      intA = product((/(z+i, i = 0, n-1)/)) / ( z + j)
    end function intA
  end function A


  function B(n, j)
    implicit none
    integer :: n, j, i
    real(8) :: B, integral
    call gauss_integral(0d0, 1d0, intB, integral)

    B = (-1d0)**(j+1) / &
        product((/(i, i=1,j+1)/)) / product((/(i, i=1,n-2-j)/)) * integral
    contains
    function intB(z)
      implicit none
      real(8) :: z, intB
      integer :: i
        intB = product((/(z + i, i = -1, n-2)/)) / ( z + j)
    end function intB
  end function B


  subroutine gauss_integral(a0, b0, func, res) ! fix n - убираем из аргументов
    implicit none
    real(8) :: a0, b0, res, func
    integer, parameter :: n = 6
    integer :: i
    real(8), dimension(1:n) :: A, t
    character(2) :: num
    interface
      function func(x)
      real(8) :: x
      real(8) :: f
      end function
    end interface

    ! Возьмём уже известные числа для n = 6
    !  n = 6

    A(1) = 0.17132449238
    t(1) = 0.93246951421
    A(2) =  A(1)
    t(2) = -t(1)
    A(3) = 0.36076157305
    t(3) = 0.66120938647
    A(4) =  A(3)
    t(4) = -t(3)
    A(5) = 0.46791399346
    t(5) = 0.23861918608
    A(6) =  A(5)
    t(6) = -t(5)

    res = 0

    do i = 1, n
        res = res + A(i) * (b0 - a0) / 2 * func(t(i)*(b0 - a0)/2 + (a0+b0)/2)
    end do

  end subroutine gauss_integral

end module coeff_adams_mod
