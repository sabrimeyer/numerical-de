function [x,y] = odeBDF2(f, fy, a, b, ya, n, sigma, tol, maxiter)
    h = (b-a)/n;
    x=a:h:b;
    m = length(ya);
    y = zeros(m, n+1);
    
    [~,y2] = thetaSchema(f, fy, a, a + h, ya, 1, 1/2, sigma, tol, maxiter);
    
    y(:,1:2) = y2(:,1:2);
    
    for i = 1:n-1
        g = @(z) z + 1/3*y(:,i) - 4/3*y(:,i+1) - 2/3*h*f(x(i+2),z);
        gp = @(z) eye(m) - 2/3*h*fy(x(i+2),z);
        y(:,i+2) = newtonIterationSWS(g, gp, y(:,i+1), sigma, tol, maxiter); %solve non-linear System
    end
    
end

