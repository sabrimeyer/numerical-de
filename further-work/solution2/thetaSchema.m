function [x, y] = thetaSchema(f, fy, a, b, ya, n, theta, sigma, tol, maxiter)
    m = length(ya);
    h = (b - a)/n;
    x = a:h:b; % n+1 Stuetzstellen
    y = zeros(m, n+1); % m Funktionen, n+1 Werte
    y(:,1) = ya; %Anfangswert
    z0 = ya;
    
    for i = 1:n
        g = @(z) z - y(:,i) - h*((1 - theta)*f(x(i), y(:,i)) + theta*f(x(i+1),z));
        gp = @(z) eye(m) - h*theta*fy(x(i+1),z);
        y(:,i+1) = newtonIterationSWS(g, gp, z0, sigma, tol, maxiter); %GS loesen
        z0 = y(:,i); %Startwertupdate
    end
    
end

