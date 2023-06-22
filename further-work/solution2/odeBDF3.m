function [x,y] = odeBDF3(f, fy, a, b, ya, n, sigma, tol, maxiter)
    h = (b-a)/n;
    x=a:h:b;
    m = length(ya);
    y = zeros(m, n+1);
    
    [~,y2] = thetaSchema(f, fy, a, a + 2*h, ya, 2, 1/2, sigma, tol, maxiter);
    
    y(:,1:3) = y2(:,1:3);
    
    for i = 1:n-2
        g = @(z) z - 2/11*y(:,i) + 9/11*y(:,i+1) - 18/11*y(:,i+2) - 6/11*h*f(x(i+3),z);
        gp = @(z) eye(m) - 6/11*h*fy(x(i+3),z);
        y(:,i+3) = newtonIterationSWS(g, gp, y(:,i+2), sigma, tol, maxiter); %solve non-linear System
    end
    
end

