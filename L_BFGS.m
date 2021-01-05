function [x,obj_val,iter] = L_BFGS(x0,eps,obj_func,gradient_func)
%% Implementation of limited memory version of BFGS Method (L-BFGS)
% 
% quasi-Newton direction:d_k = -H_k*gradient(x_k)
% two-loop recursion of L-BFGS (approximated H_k*gradient(x_k))
%   - s_i       = x_{i+1}-x_i  ;  
%   - y_i       = gradient(x_{i+1})-gradient(x_i))
%   - H0_k      = gamma_k*I      
%   - gamma_k   = ((s_{k-1})^T*y_{k-1})/((y_{k-1})^T*y_{k-1}) 
%   - (or)      = ((s_{k-1})^T*s_{k-1})/((s_{k-1})^T*y_{k-1})
%   - m         : memory
%   Set q = gradient(x_k)
%   for i = k-1,k-2,...,k-m do
%       rho_i = 1/((s_k)^T*y_k)
%       a_i   = rho_i*(s_i)^T*q
%       q     = q-a_i*y_i
%   Set r = H0_k*q
%   for i = k-m,k-m+1,...,k-1 do
%       beta  = rho_i*(y_i)^T*r
%       r     = r + (a_i-beta)*beta_i
%   STOP with result r = H_k*gradient(x_k)
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

%% parameter for L-BFGS
m = 5;

%% buffer
n = length(x0);
sk_struct = struct("list",zeros(n,m),"num",0);  % store s_k
yk_struct = struct("list",zeros(n,m),"num",0);  % store y_k

%% iterative process
while ng > eps && iter < MAX_ITER
    iter    = iter + 1;

    sk_list = sk_struct.list;
    yk_list = yk_struct.list;
    cur_m   = sk_struct.num;
    d_k     = -g;
    if iter ~= 1
        d_k = -L_BFGS_loop(g,sk_list,yk_list,cur_m);
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

    % update memoruy buffer
    sk    = x_new - x_k;
    g_new = gradient_func(x_new);
    yk    = g_new - g;
    if cur_m < m
        sk_struct.list(:,cur_m+1) = sk;
        yk_struct.list(:,cur_m+1) = yk;
        sk_struct.num = cur_m+1;
        yk_struct.num = cur_m+1;
    else
        sk_struct.list(:,1:m-1) = sk_struct.list(:,2:m);
        sk_struct.list(:,m) = sk;
        yk_struct.list(:,1:m-1) = yk_struct.list(:,2:m);
        yk_struct.list(:,m) = yk;
    end
    
    % update
    x_k     = x_new;
    obj_val = obj_func(x_k);
    g       = g_new;
    ng      = norm(g);
end
x = x_k;

%% Implemetation of two-loop recursion in L-BFGS
function r = L_BFGS_loop(g,sk_list,yk_list,m)
    function A = identity( N )
        A(N,N) = 0;
        A((N+1)*(0:N-1)+1) = 1;
    end
    q = g;
    alpha = zeros(m,1);
    rho = zeros(m,1);
    for i = m:1
        si = sk_list(:,i);
        yi = yk_list(:,i);
        rho(i) = 1/(transpose(si)*yi);
        alpha(i) = rho(i)*transpose(si)*q;
        q = q-alpha(i)*yi;
    end
    I = identity(length(g));
    gamma_k = (transpose(sk_list(:,m))*yk_list(:,m))/(transpose(yk_list(:,m))*yk_list(:,m));
    H_k_0 = gamma_k*I;
    r = H_k_0 * q;
    for i = 1:m
        si = sk_list(:,i);
        yi = yk_list(:,i);
        beta = rho(i)*transpose(yi)*r;
        r = r + (alpha(i)-beta)*si;
    end
end
end

