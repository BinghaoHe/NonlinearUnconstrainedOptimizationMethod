function V = matrix2vector(M)
[m,n] = size(M);
V = zeros(m+n,1);
for i = 1:m
    for j = 1:n
        V((i-1)*n+j) = M(i,j);
    end
end

end

