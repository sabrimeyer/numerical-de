function z = newtonIterationSWS(g, gp, z0, sigma, tol, maxiter)
    z = z0;
    
    for i = 1:maxiter
        d = gp(z)\g(z); % GS l?sen
        if norm(d) <= tol
            break % Nullstelle gefunden
        end
        alpha = 1;
        while norm(g(z - alpha.*d)) > norm(g(z))*(1 - sigma*alpha)
            alpha = alpha/2;
        end
        z = z - alpha*d;
    end
    
end

