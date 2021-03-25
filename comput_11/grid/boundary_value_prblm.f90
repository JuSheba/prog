program boundary_value_prblm
use grid_method_mod
use param_func
implicit none
  real(8) :: res
    res  = grid(func, a, b, alpha, beta, num_int)

end program boundary_value_prblm
