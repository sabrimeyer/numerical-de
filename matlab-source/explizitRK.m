function [t, y] = explizitRK(F, t0, y0, T, n, methode)
% Initialisierung
    h = (T-t0)./n; %Schrittweite
    m = length(y0); %Dimension
    t = zeros(n+1,1); %Zeitintervall
    t(1) = t0;
    y = zeros(m,n); % m Dimensionen, n Werte
    y(:,1) = y0;
    
    
    switch methode %Fallunterscheidung der Verfahren
        case 'Euler'
            for i = 1:n
                t(i+1) = t(i)+h;
                y(:,i+1) = y(:,i) + h*F(t(i),y(:,i));
            end
            
        case 'Heun2'
            for i = 1:n
                t(i+1) = t(i)+h;
                k1 = F(t(i),y(:,i));
                k2 = F(t(i)+h,y(:,i)+h*k1);
                y(:,i+1) = y(:,i) + h/2*(k1+k2);
            end
            
        case 'Kutta3'
            for i = 1:n
                t(i+1) = t(i)+h;
                k1 = F(t(i),y(:,i));
                k2 = F(t(i)+h/2,y(:,i)+h/2*k1);
                k3 = F(t(i)+h,y(:,i)-h*k1+2*h*k2);
                y(:,i+1) = y(:,i) + h/6*(k1+4*k2+k3);
            end
        case 'RK4'
            for i = 1:n
                k1 = F(t(i),y(:,i));
                k2 = F(t(i)+h/2,y(:,i)+h/2*k1);
                k3 = F(t(i)+h/2,y(:,i)+h/2*k2);
                k4 = F(t(i)+h,y(:,i)+h*k3);
                t(i+1) = t(i)+h;
                y(:,i+1) = y(:,i) + h/6*(k1+2*k2+2*k3+k4);
            end
    end
    
    
end

