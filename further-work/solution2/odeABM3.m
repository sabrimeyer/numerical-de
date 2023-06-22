function [x,y] = odeABM3(f, a, b, ya, n, m)
h = (b - a)/n;
x = a:h:b;
m = length(ya);
y = zeros(m, n+1);

% Startwerte berechnen
[~, y2] = klassischRK(f, a, ya, h);
y(:,1:3) = y2(:,1:3);

beta = [5/12, -4/3, 23/12, 0];
yy = y;
fval = [f(x(1),y(:,1)), f(x(2),y(:,2)), f(x(3),y(:,3))];

for i = 1:n-2
    yy(:,i+3) = y(:,i+2) + h*beta(1)*fval(:,1) + h*beta(2)*fval(:,2) + h*beta(3)*fval(:,3); % Prädiktor-Schritt
    for j = 1:m
        temp = f(x(i+3),yy(:,i+3));
        yy(:,i+3) = y(:,i+2) + h*1/24*fval(:,1) - h*5/24*fval(:,2) + h*19/24*fval(:,3) + h*3/8*temp; % Korrektor-Schritt
    end
    % Update
    fval = [fval(:,2), fval(:,3), temp];
    y(:,i+3) = yy(:,i+3);
end
y = yy;

end

function [x2, y2] = klassischRK(f, a, ya, h)
    alpha2 = [0, 1/2, 1/2, 1];
    beta2 = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
    gamma2 = [1/6 1/3 1/3 1/6];
    [x2, y2] = explizitRK(f, a, a + 2*h, ya, 2, alpha2, beta2, gamma2);
end


