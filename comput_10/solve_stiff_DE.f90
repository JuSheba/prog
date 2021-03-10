program solve_stiff_DE
use param_func_mod
use stiff_DE_mod
implicit none

real(8):: result(size(X0))
  result = Rosenbrock(func, x0, h, t_0, t_k)
  result = PC_Method(func, x0, h, t_0, t_k)

end program solve_stiff_DE
