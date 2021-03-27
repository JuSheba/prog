module grid_method_mod
  use newton_mod
  use param_func
  implicit none

contains

  function grid(func, a, b, alpha, beta, num_int)
    implicit none
    interface
      function func(t, x, dx)
        implicit none
        real(8), dimension(:) :: x, dx
        real(8), dimension(size(x)) :: func
        real(8) :: t
      end function func
    end interface
    real(8), dimension(1:) :: alpha, beta
    real(8), dimension(1:size(alpha)) :: x_0, x_k
    integer :: num_int, n, i, j
    real(8) :: t_0, t_k, h1, a, b, grid, t
    real(8), dimension(1:num_int) :: time!, solution
    real(8), dimension(1:size(alpha)*(num_int + 1)) :: all_x
    real(8), dimension(1:size(alpha)*(num_int - 1)) :: x_res, zero

      n = size(alpha)
      h1 = (b - a) / num_int
      x_0 = alpha; x_k = beta
      t_0 = a;     t_k = b

      zero = 0d0
      x_res(1:n*(num_int-1)) = 0d0
      forall(i = 1:num_int) time(i) = a + i * h1

      call newton(zero, 1000, x_res, get_grid)

      ! это если хочется сверить с решением
      !do i = 1, num_int
      !   solution(i) = check_result(time(i))
      !end do
      open(22, file = 'grid.dat')
      write(22,*) a, x_0
        do i = 1, num_int - 1
          write(22,*) time(i), x_res(n*(i-1)+1:n*i)!, solution(i)
        end do
      write(22,*) b, x_k
      close(22)

contains

    function get_grid(x_res) result(res)
      real(8), dimension(:) :: x_res, res(1:size(x_res))
      real(8), dimension(1:size(x_0)) :: x1, x2, x3, dx
      real(8), dimension(1:size(x_0)*(num_int+1)) :: all_x
        all_x(1:n) = x_0
        all_x(size(all_x)-n+1 : size(all_x)) = x_k
        all_x(n+1 : size(all_x)-n) = x_res

        do i = 1, num_int - 1
            x1 = all_x(n*(i-1)+1 : n*i)
            x2 = all_x(n*i+1     : n*(i+1))
            x3 = all_x(n*(i+1)+1 : n*(i+2))
            dx = func(time(i), x2, (x3 - x1) / (2*h1))
            do j = 1, n
                res((i-1)*n + j) = x1(j) - 2*x2(j) + x3(j) - h1**2 * dx(j)
    		end do
        end do
    end function get_grid

  end function grid

end module grid_method_mod
