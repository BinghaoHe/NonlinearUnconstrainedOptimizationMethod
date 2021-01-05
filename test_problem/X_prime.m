function X_prime = X_prime(X,boundary_f)
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

end

