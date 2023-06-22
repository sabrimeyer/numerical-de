function [x, y] = explizitRKadapt(f, a, b, ya, h0, alpha, beta, gamma, gammah, p, tau, epsilon)
    m = length(alpha);
    N = length(ya);
    h = h0;
    y(:,1) = ya; %Anfangswert
    k = zeros(N, m);
    
    x = a;
    i = 1;
    delta = epsilon;

    while x(i) < b
        k(:,1) = f(x(i), y(:,i));
        ySum = gamma(1)*k(:,1);
        deltaSum = (gamma(1) - gammah(1))*k(:,1);
        
        for j = 2:m %outer sum
            kSum = 0;
            for l = 1:j-1 %inner sum
                kSum = kSum + beta(j,l)*k(:,l);
            end
            k(:,j) = f(x(i) + h*alpha(j), y(:,i) + h*kSum);
            ySum = ySum + gamma(j)*k(:,j);
            deltaSum = (gamma(j) - gammah(j))*k(:,j);
        end
        
        y = [y, zeros(N,1)];
        while 1
            h = tau*(epsilon/delta)^(1/(p+1))*h;
            y(:,i+1) = y(:,i) + h*ySum;
            delta = h*norm(deltaSum);
            
            if delta <= epsilon
                break
            end
        end
        
            x = [x, x(i) + h]; %x(i+1) = x(i) + h;
            i = i + 1;
    end

end

