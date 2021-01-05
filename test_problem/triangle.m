function area = triangle(X,r)
    x1 = X(1);
    x2 = X(2);
    x3 = X(3);
    
    area = (r/2)*sqrt((x1-x2)^2+(x1-x3)^2+r^2);
end

