function dydt = geodesic(t,y)
    % This script encode the right hand side of the ODE as a function.
    % Using f(x_1,x_2) = 1 + x_1
    % Need to set alpha manually
    % Note that y_1 = x_1, y_2 = x_2, y_3 = dx_1/dt and y_4 = dx_2/dt
    
    p = 3; d = 2;
    alpha = -2*(p-1)/d;
    f_inv = (1+y(1))^(-1);
    grad_f = [1; 0];
    
    dydt = zeros(4,1);
    dydt(1) = y(3);
    dydt(2) = y(4);
    dydt(3) = -alpha*f_inv*y(3)^2 + (alpha/2)*f_inv*(y(3)^2 + y(4)^2);
    dydt(4) = -alpha*f_inv*y(3)*y(4);
    
end

    
    
