function [x,y] = odeABM2(f, a, b, ya, n, m)
h = (b - a)/n;
x = a:h:b;
m = length(ya);
y = zeros(m, n+1);

% Startwerte berechnen
[~, y2] = klassischRK(f, a, ya, h);
y(:,1:2) = y2(:,1:2);

beta = [-1/2, 3/2, 0];
yy = y;
fval = [f(x(1),y(:,1)), f(x(2),y(:,2))];

for i = 1:n-1
    yy(:,i+2) = y(:,i+1) + h*beta(1)*fval(:,1) + h*beta(2)*fval(:,2); % Prädiktor-Schritt
    for j = 1:m
        temp = f(x(i+2),yy(:,i+2));
        yy(:,i+2) = y(:,i+1) - h*1/12*fval(:,1) + h*2/3*fval(:,2) + h*5/12*temp; % Korrektor-Schritt
    end
    % Update
    fval = [fval(:,2), temp];
    y(:,i+2) = yy(:,i+2);
end
y = yy;

end

function [x2, y2] = klassischRK(f, a, ya, h)
    alpha2 = [0, 1/2, 1/2, 1];
    beta2 = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
    gamma2 = [1/6 1/3 1/3 1/6];
    [x2, y2] = explizitRK(f, a, a + h, ya, 1, alpha2, beta2, gamma2);
end
