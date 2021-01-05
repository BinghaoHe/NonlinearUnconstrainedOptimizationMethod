%% This is a testing script for different optimization problem
% The tested problem is the discretized minimum surface problem
addpath(genpath(pwd));
obj_func        = @obj; 
gradient_func   = @obj_gradient; 
hessian_func    = @obj_hessian; 
k               = 21; % degree of discretization in minimum surface problem

%% Setting up
x0  = rand((k-2)*(k-2),1);      % random initialization
eps = 1e-6;                     % set the 

%% Invoke different optimization method %%
[x,opt_val,iter] = backtracking(x0,eps,obj_func,gradient_func);
[x,opt_val,iter] = newton(x0,eps,obj_func,gradient_func);
[x,opt_val,iter] = L_BFGS(x0,eps,obj_func,gradient_func);

%% Plot the result of the tested minimum surface problem
figure;
A       = vector2matrix(x);
A_prime = X_prime(A,@f1);
a       = matrix2vector(A_prime);
k       = sqrt(length(a));
z       = zeros(k,k);
for i = 1:k
    z(i,:)=a((i-1)*k+1:i*k);
end
x = 0:1/(k-1):1;
y = 0:1/(k-1):1;
[x,y] = meshgrid(x,y);
T = delaunay(x,y);
trisurf(T,x,y,z);
title("Minimum Surface");
xlabel("x");
ylabel("y");

set(gcf,'position',[100,100,500,500])

