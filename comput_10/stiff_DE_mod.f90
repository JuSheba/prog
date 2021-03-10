module stiff_DE_mod
  use elimination_mod
  use param_func_mod
  use methods_adams_mod
  implicit none
  contains
!__________________________________________________________________________________
  function Rosenbrock(func, x0, h, t_0, t_k)
    implicit none
    real(8) :: alpha, beta, gamma, eps, t, t_0, t_k, h
    real(8), dimension(:)   :: x0
    real(8), dimension(1:size(X0)) :: x, buf, system_vctr, Rosenbrock, res
    real(8), dimension(1:size(X0),1:size(X0)) :: jacobian, E, system_mtrx
    integer :: i, j
    interface
      function func(t, x)
        real(8), dimension(:) :: X
        real(8), dimension(1:size(X0)) :: func
        real(8) :: t
      end function func
    end interface
    !_______________________________________________________________________
    alpha =  1.077d0
    beta  = -0.372d0
    gamma = -0.577d0
    eps = epsilon(1d0)
    t = t_0

    E = 0d0
    forall(i = 1:size(x0)) E(i,i) = 1d0
    Rosenbrock = x0

    open(11, file = 'rosen.dat')
    do while (t < t_k)
      buf = Rosenbrock
      do i = 1, size(x0)
          buf = Rosenbrock
          buf(i) = Rosenbrock(i) + eps
          jacobian(:,i) = (func(t, buf) - func(t, Rosenbrock)) / eps
      end do
      !write(6,*) jacobian

      system_mtrx = E - alpha * h * jacobian - beta * h**2 * matmul(jacobian,jacobian)
      system_vctr = func(t, gamma * h * func(t, Rosenbrock) + Rosenbrock) * h
      call Gauss(system_mtrx, system_vctr, res)
      Rosenbrock = Rosenbrock + res
      write(11,*) t, Rosenbrock
      t = t + h
    end do
    close(11)
  end function Rosenbrock
!__________________________________________________________________________________
  function PC_Method(func, x0, h, t_0, t_k)
    implicit none
    real(8) :: t, t_0, t_k, h
    real(8), dimension(:)   :: x0
    real(8), dimension(1:size(X0)) :: PC_Method, res, x
    integer :: i, j
    interface
      function func(t, x)
        real(8), dimension(:) :: x
        real(8), dimension(1:size(x)) :: func
        real(8) :: t
      end function func
    end interface
!    !_______________________________________________________________________
    t = t_0
    PC_Method = x0

    open(22, file = 'precor.dat')
      do while (t < t_k)
        res = e_adams(PC_Method, t, x)
        PC_Method = i_adams(res, t, x)
        t = t + h
        write(22,*) t, PC_Method
      end do
    close(22)
  end function PC_Method


end module stiff_DE_mod
