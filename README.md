# Nonlinear Unconstrained Optimization Methods
This repository contrains three different methods implemented with MatLab for nonlinear unconstrained optimization, including Gradient Descent with Backtracking, Globalized Newton, and L-BFGS methods.


## Gradient descent method with Armijo line search (backtracking) step size strategy
backtracking.m (function)
Input:
  x0:             - the initial value of x
  eps:            - the tolerance for norm of the gradient
                      If norm(gradient) < eps, then STOP
  obj_func        - a pointer to objective function
                      obj_func(x): scalar value
  gradient_func   - a pointer to gradient function
                      gradient(x): a vector
Output:
  x:              - the solution point
  obj_val         - the objective value under the solution point
  iter            - the number of iterations carried
  
## Globalized Newto method
newton.m (function)
Input:
  x0:             - the initial value of x
  eps:            - the tolerance for norm of the gradient
                      If norm(gradient) < eps, then STOP
  obj_func        - a pointer to objective function
                      obj_func(x): scalar value
  gradient_func   - a pointer to gradient function
                      gradient(x): a vector
  hessian_func    - a pointer to hessian function
                      gradient(x): a matrix
Output:
  x:              - the solution point
  obj_val         - the objective value under the solution point
  iter            - the number of iterations carried
  
## Limited memory BFGS (L-BFGS) method
L_BFGS.m (function)
Input:
  x0:             - the initial value of x
  eps:            - the tolerance for norm of the gradient
                      If norm(gradient) < eps, then STOP
  obj_func        - a pointer to objective function
                      obj_func(x): scalar value
  gradient_func   - a pointer to gradient function
                      gradient(x): a vector
Output:
  x:              - the solution point
  obj_val         - the objective value under the solution point
  iter            - the number of iterations carried

## Example
test.m (script)
Use the discretized minimum surface problem as an example to show how to use different optimization methods.
