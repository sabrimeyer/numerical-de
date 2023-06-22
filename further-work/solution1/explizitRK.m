function [x, y] = explizitRK(f, a, b, ya, n, alpha, beta, gamma)
    m = length(alpha);
    N = length(ya);
    h = (b - a)/n; %Schrittweite
    x = a:h:b; %Stuetzstellen
    y = zeros(N, n+1);
    y(:,1) = ya; %Anfangswert
    k = zeros(N, m);
    
    for i = 1:n
        k(:,1) = f(x(i), y(:,i));
        ySum = gamma(1)*k(:,1);
        
        for j = 2:m %outer sum
            kSum = 0;
            for l = 1:j-1 %inner sum
                kSum = kSum + beta(j,l)*k(:,l);
            end
            k(:,j) = f(x(i) + h*alpha(j), y(:,i) + h*kSum);
            ySum = ySum + gamma(j)*k(:,j);
        end
        
        y(:,i+1) = y(:,i) + h*ySum;
    end
    
end

