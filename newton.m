function [x,obj_val,iter] = newton(x0,eps,obj_func,gradient_func,hessian_func)
%% Implementation of Globalized Newton Method 
% with Armijo line search (backtracking) step size strategy
% 
% Newton direction: s_k  s.t. hessian(x_k)*s_k = -gradient(x_k)
% If Newton direction dosen't exists,
% or -gradient(x_k)^T*s_k < min{beta1,beta2*norm(s_k)}*norm(s_k)^2
%   descent direction d_k = -gradient(x_k)
% Otherwise, d_k = s_k
%
% Armijo line search (backtracking) is used to calculate the step size
%
% Input:
%   x0:             - the initial value of x
%   eps:            - the tolerance for norm of the gradient
%                       If norm(gradient) < eps, then STOP
%   obj_func        - a pointer to objective function
%                       obj_func(x): scalar value
%   gradient_func   - a pointer to gradient function
%                       gradient(x): a vector
%   hessian_func    - a pointer to hessian function
%                       gradient(x): a matrix
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
h         = hessian_func(x_k);            % hessian(x_k)
obj_val   = obj_func(x_k);                % objective value(x_k)
ng        = norm(g);                      % norm of gradient  

%% parameters for globalized newton method
beta1 = 1e-6;
beta2 = 1e-6;

%% parameters for Armijo line search (backtracking) step size strategy
alpha_0   = 1;                            % the initial step size
sigma     = 0.5;
gamma     = 1e-2; 

%% iterative process
while ng > eps && iter < MAX_ITER
    iter  = iter + 1;
    d_k = -g;
    if det(h) ~= 0          % hessian is invertible
        s_k = -(h^-1)*g;    % Newton direction
        if -transpose(g)*s_k >= min(beta1,beta2*norm(s_k))*(norm(s_k)^2)
            d_k = s_k;
        end
    end
    
    % step size strategy
    alpha   = alpha_0;
    x_new   = x_k + alpha*d_k;
    obj_new = obj_func(x_new);
    r       = -gamma*ng^2;
    % find alpha satisfying Armijo condition
    while obj_new - obj_val > r*alpha
        alpha       = sigma*alpha;
        x_new       = x_k + alpha*d_k;
        obj_new     = obj_func(x_new);
    end

    % update
    x_k     = x_new;
    obj_val = obj_func(x_k);
    g       = gradient_func(x_k);
    h       = hessian_func(x_k);
    ng      = norm(g);
end
x = x_k ;
end

