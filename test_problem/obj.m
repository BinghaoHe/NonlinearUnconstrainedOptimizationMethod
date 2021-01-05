function obj_val = obj(x)
boundary_f = @f1;
X = vector2matrix(x);
%% X (k-2)*(k-2)
size_X = size(X);
k = size_X(1)+2; %% k*k, the grid is (k-1)*(k-1)
r = 1/(k-1);

%% Form new X_prime k*k
X_prime = zeros(k,k);
for i = 1:k
    for j = 1:k
        if i == 1 || i == k || j == 1 || j == k %% boundary  x: |   y:->
            temp_x = (i-1)*(1/(k-1));
            temp_y = (j-1)*(1/(k-1));
            X_prime(i,j) = boundary_f(temp_x,temp_y);
        else
            X_prime(i,j) = X(i-1,j-1);
        end
    end
end

%% Calculate the triangles' area
obj_val = 0;
for i = 1:k-1
    for j = 1:k-1
       %% up triangle
        x1 = X_prime(i,j+1);
        x2 = X_prime(i,j);
        x3 = X_prime(i+1,j+1);
        obj_val = obj_val + triangle([x1,x2,x3],r);
       %% down triangle
        x1 = X_prime(i+1,j);
        x2 = X_prime(i,j);
        x3 = X_prime(i+1,j+1);
        obj_val = obj_val + triangle([x1,x2,x3],r);    
    end
end

end

