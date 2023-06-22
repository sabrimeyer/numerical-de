function [x, y] = odeAB3(f, a, b, ya, n)
    
    h = (b - a)/n;
    x = a:h:b;
    m = length(ya);
    y = zeros(m, n+1);
    
    beta = [5/12, -4/3, 23/12, 0];
    [~, y2] = klassischRK(f, a, b, ya, h);
    y(:,1:3) = y2(:,1:3);
    
    fval = [f(x(1),y(:,1)), f(x(2),y(:,2)), f(x(3), y(:,3))]; % k=3
    
    for i = 1:n-2
        y(:,i+3) = y(:,i+2) + h*beta(1)*fval(:,1) + h*beta(2)*fval(:,2) + h*beta(3)*fval(:,3);
        fval = [fval(:,2), fval(:,3), f(x(i+3),y(:,i+3))];
    end

end

function [x2, y2] = klassischRK(f, a, b, ya, h)
    alpha2 = [0, 1/2, 1/2, 1];
    beta2 = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
    gamma2 = [1/6 1/3 1/3 1/6];
    [x2, y2] = explizitRK(f, a, a + 2*h, ya, 2, alpha2, beta2, gamma2);
end



