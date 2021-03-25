program boundary_value_prblm
use shooting_method_mod
use param_func
implicit none

  real(8) :: res
    res = shooting(func, a, b, alpha, beta, num_int, 'left')

end program boundary_value_prblm
