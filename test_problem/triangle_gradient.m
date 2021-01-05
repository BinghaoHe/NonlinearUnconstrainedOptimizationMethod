function gradient = triangle_gradient(X,r)
    x1 = X(1);
    x2 = X(2);
    x3 = X(3);
    
    df_x1 = -(r*(2*x2 - 4*x1 + 2*x3))/(4*((x1 - x2)^2 + (x1 - x3)^2 + r^2)^(1/2));
    df_x2 = -(r*(2*x1 - 2*x2))/(4*((x1 - x2)^2 + (x1 - x3)^2 + r^2)^(1/2));
    df_x3 = -(r*(2*x1 - 2*x3))/(4*((x1 - x2)^2 + (x1 - x3)^2 + r^2)^(1/2));
    gradient = [df_x1;df_x2;df_x3];
end

