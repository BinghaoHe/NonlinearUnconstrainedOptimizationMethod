function hessian = obj_hessian(x)
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

%% Calculate the gradient
hessian = zeros((k-2)*(k-2),(k-2)*(k-2));

    %% up triangle
    %  2 ！！！！！！   1
    %       -       |
    %            -  |
    %               | 3
for i = 1:k-1
    for j = 1:k-1
        x1 = X_prime(i,j+1);
        x2 = X_prime(i,j);
        x3 = X_prime(i+1,j+1);
        temp_hessian = triangle_hessian([x1,x2,x3],r);
        if i == 1 || i == k-1 || j == 1 || j == k-1
            if i == 1 && j == 1      %% 4 corners
                hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j) = hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j) + temp_hessian(3,3);
            elseif i == 1 && j == k-1
                0;
            elseif i == k-1 && j == 1
                hessian((i-2)*(k-2)+j,(i-2)*(k-2)+j) = hessian((i-2)*(k-2)+j,(i-2)*(k-2)+j) + temp_hessian(1,1);
            elseif i == k-1 && j == k-1
                hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j-1) = hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j-1) + temp_hessian(2,2);
            elseif i == 1  %% 4 edges
                hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j) = hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j) + temp_hessian(3,3);
            elseif i == k-1
                hessian((i-2)*(k-2)+j,(i-2)*(k-2)+j) = hessian((i-2)*(k-2)+j,(i-2)*(k-2)+j) + temp_hessian(1,1);
                hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j-1) = hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j-1) + temp_hessian(2,2);
                hessian((i-2)*(k-2)+j,(i-2)*(k-2)+j-1) = hessian((i-2)*(k-2)+j,(i-2)*(k-2)+j-1) + temp_hessian(1,2);
                hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j) = hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j) + temp_hessian(2,1);
            elseif j == 1            
                hessian((i-2)*(k-2)+j,(i-2)*(k-2)+j) = hessian((i-2)*(k-2)+j,(i-2)*(k-2)+j) + temp_hessian(1,1);
                hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j) = hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j) + temp_hessian(3,3);
                hessian((i-2)*(k-2)+j,(i-1)*(k-2)+j) = hessian((i-2)*(k-2)+j,(i-1)*(k-2)+j) + temp_hessian(1,3);
                hessian((i-1)*(k-2)+j,(i-2)*(k-2)+j) = hessian((i-1)*(k-2)+j,(i-2)*(k-2)+j) + temp_hessian(3,1);
            elseif j == k-1   
                hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j-1) = hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j-1) + temp_hessian(2,2);
            end
        else         
            hessian((i-2)*(k-2)+j,(i-2)*(k-2)+j) = hessian((i-2)*(k-2)+j,(i-2)*(k-2)+j) + temp_hessian(1,1);
            hessian((i-2)*(k-2)+j,(i-2)*(k-2)+j-1) = hessian((i-2)*(k-2)+j,(i-2)*(k-2)+j-1) + temp_hessian(1,2);
            hessian((i-2)*(k-2)+j,(i-1)*(k-2)+j) = hessian((i-2)*(k-2)+j,(i-1)*(k-2)+j) + temp_hessian(1,3);
            hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j) = hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j) + temp_hessian(2,1);
            hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j-1) = hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j-1) + temp_hessian(2,2);
            hessian((i-2)*(k-2)+j-1,(i-1)*(k-2)+j) = hessian((i-2)*(k-2)+j-1,(i-1)*(k-2)+j) + temp_hessian(2,3);
            hessian((i-1)*(k-2)+j,(i-2)*(k-2)+j) = hessian((i-1)*(k-2)+j,(i-2)*(k-2)+j) + temp_hessian(3,1);
            hessian((i-1)*(k-2)+j,(i-2)*(k-2)+j-1) = hessian((i-1)*(k-2)+j,(i-2)*(k-2)+j-1) + temp_hessian(3,2);
            hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j) = hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j) + temp_hessian(3,3);
       end
    end
end

%% down triangle
    %  2 
    %     |  -       
    %     |       -  
    %  1  |！！！！！！  3
for i = 1:k-1
    for j = 1:k-1
        x1 = X_prime(i+1,j);
        x2 = X_prime(i,j);
        x3 = X_prime(i+1,j+1);
        temp_hessian = triangle_hessian([x1,x2,x3],r);
        if i == 1 || i == k-1 || j == 1 || j == k-1
            if i == 1 && j == 1      %% 4 corners
                hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j) = hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j) + temp_hessian(3,3);
            elseif i == 1 && j == k-1
                hessian((i-1)*(k-2)+j-1,(i-1)*(k-2)+j-1) = hessian((i-1)*(k-2)+j-1,(i-1)*(k-2)+j-1) + temp_hessian(1,1);
            elseif i == k-1 && j == 1
                0;
            elseif i == k-1 && j == k-1
                hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j-1) = hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j-1) + temp_hessian(2,2);
            elseif i == 1  %% 4 edges
                hessian((i-1)*(k-2)+j-1,(i-1)*(k-2)+j-1) = hessian((i-1)*(k-2)+j-1,(i-1)*(k-2)+j-1) + temp_hessian(1,1);
                hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j) = hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j) + temp_hessian(3,3);
                hessian((i-1)*(k-2)+j-1,(i-1)*(k-2)+j) = hessian((i-1)*(k-2)+j-1,(i-1)*(k-2)+j) + temp_hessian(1,3);
                hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j-1) = hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j-1) + temp_hessian(3,1);
            elseif i == k-1
                hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j-1) = hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j-1) + temp_hessian(2,2);
            elseif j == 1
                hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j) = hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j) + temp_hessian(3,3);
            elseif j == k-1
                hessian((i-1)*(k-2)+j-1,(i-1)*(k-2)+j-1) = hessian((i-1)*(k-2)+j-1,(i-1)*(k-2)+j-1) + temp_hessian(1,1);
                hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j-1) = hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j-1) + temp_hessian(2,2);
                hessian((i-1)*(k-2)+j-1,(i-2)*(k-2)+j-1) = hessian((i-1)*(k-2)+j-1,(i-2)*(k-2)+j-1) + temp_hessian(1,2);
                hessian((i-2)*(k-2)+j-1,(i-1)*(k-2)+j-1) = hessian((i-2)*(k-2)+j-1,(i-1)*(k-2)+j-1) + temp_hessian(2,1);
            end
        else
            hessian((i-1)*(k-2)+j-1,(i-1)*(k-2)+j-1) = hessian((i-1)*(k-2)+j-1,(i-1)*(k-2)+j-1) + temp_hessian(1,1);
            hessian((i-1)*(k-2)+j-1,(i-2)*(k-2)+j-1) = hessian((i-1)*(k-2)+j-1,(i-2)*(k-2)+j-1) + temp_hessian(1,2);
            hessian((i-1)*(k-2)+j-1,(i-1)*(k-2)+j) = hessian((i-1)*(k-2)+j-1,(i-1)*(k-2)+j) + temp_hessian(1,3);
            hessian((i-2)*(k-2)+j-1,(i-1)*(k-2)+j-1) = hessian((i-2)*(k-2)+j-1,(i-1)*(k-2)+j-1) + temp_hessian(2,1);
            hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j-1) = hessian((i-2)*(k-2)+j-1,(i-2)*(k-2)+j-1) + temp_hessian(2,2);
            hessian((i-2)*(k-2)+j-1,(i-1)*(k-2)+j) = hessian((i-2)*(k-2)+j-1,(i-1)*(k-2)+j) + temp_hessian(2,3);
            hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j-1) = hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j-1) + temp_hessian(3,1);
            hessian((i-1)*(k-2)+j,(i-2)*(k-2)+j-1) = hessian((i-1)*(k-2)+j,(i-2)*(k-2)+j-1) + temp_hessian(3,2);
            hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j) = hessian((i-1)*(k-2)+j,(i-1)*(k-2)+j) + temp_hessian(3,3);
       end
    end
end


end

