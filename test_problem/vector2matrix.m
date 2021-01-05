function M = vector2matrix(V)
k = sqrt(length(V));
M = zeros(k,k);
for i = 1:k
    for j = 1:k
        M(i,j) = V((i-1)*k+j);
    end
end
end

