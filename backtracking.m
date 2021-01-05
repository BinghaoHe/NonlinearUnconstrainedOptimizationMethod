function [x,obj_val,iter] = backtracking(x0,eps,obj_func,gradient_func)
%% Implementation of Gradient Descent Method 
% with Armijo line search (backtracking) step size strategy
% 
% Input:
%   x0:             - the initial value of x
%   eps:            - the tolerance for norm of the gradient
%                       If norm(gradient) < eps, then STOP
%   obj_func        - a pointer to objective function
%                       obj_func(x): scalar value
%   gradient_func   - a pointer to gradient function
%                       gradient(x): a vector
%
% Output:
%   x:              - the solution point
%   obj_val         - the objective value under the solution point
%   iter            - the number of iterations carried
%
%% initialization
MAX_ITER  = 1000;                         % maximum number of iterations
iter      = 0;                            % k
x_k       = x0;                           % x_k
g         = gradient_func(x_k);           % gradient(x_k)
obj_val   = obj_func(x_k);                % objective value(x_k)
ng        = norm(g);                      % norm of gradient  

%% parameters for Armijo line search (backtracking) step size strategy
alpha_0   = 1;                            % the initial step size
sigma     = 0.5;
gamma     = 1e-2; 

%% iterative process
while ng > eps && iter < MAX_ITER
    iter    = iter + 1;

    alpha   = alpha_0;
    x_new   = x_k - alpha*g;
    obj_new = obj_func(x_new);
    r       = -gamma*ng^2;
    % find alpha satisfying Armijo condition
    while obj_new - obj_val > r*alpha
        alpha       = sigma*alpha;
        x_new       = x_k - alpha*g;
        obj_new     = obj_func(x_new);
    end

    % update
    x_k     = x_new;
    obj_val = obj_func(x_k);
    g       = gradient_func(x_k);
    ng      = norm(g);
end
x = x_k;
end

