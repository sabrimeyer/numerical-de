function [x,y] = odeAM3(f, fy, a, b, ya, n, sigma, tol, maxiter)
    h = (b-a)/n;
    x=a:h:b;
    m = length(ya);
    y = zeros(m, n+1);
    
    beta = [1/24, -5/24, 19/24, 3/8];
    [~,y2] = thetaSchema(f, fy, a, a + 2*h, ya, 2, 1/2, sigma, tol, maxiter);
    
    y(:,1:3) = y2(:,1:3);
    
    fval = [f(x(1),y(:,1)), f(x(2),y(:,2)), f(x(3),y(:,3))]; % k=3
    
    for i = 1:n-2
        g = @(z) z - y(:,i+2) - h*beta(1)*fval(:,1) - h*beta(2)*fval(:,2) - h*beta(3)*fval(:,3) - h*beta(4)*f(x(i+3),z);
        gp = @(z) eye(m) - h*beta(4)*fy(x(i+3),z);
        y(:,i+3) = newtonIterationSWS(g, gp, y(:,i+2), sigma, tol, maxiter); %solve non-linear System

        fval = [fval(:,2), fval(:,3), f(x(i+3),y(:,i+3))];
        
    end
    
end

