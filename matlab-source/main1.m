clc
clear all
close all

% Initialisierung: Intervall [t0,T], Anfangswert y0
t0 = 1;
y0 = [4 2 0].';
T = 3;

% Initisalisierung: Funktion F, analytische Lsg y_exact
F = @(t,y) [atan(y(1)^3+2*t)+cos(y(3)), y(2)^2-t*y(1), 1/t*exp(atan(y(3)^2))+exp(-y(2)^2)+13*cos(4*pi*t)].';
y_exact = [6.986845270192330; -4.461339962786163; 3.28773579037];
err = zeros(12,4); %
h = zeros(1,12); %

for i = 2:13
    n = 2^i; %2^2, 2^3, ..., 2^13
    h(i-1) = (T-t0)/n; % Schrittweite
    
    % Loesungsverfahren [Zeitpunkt t, Approximation y]
    [t,yEuler] = explizitRK(F, t0, y0, T, n, 'Euler');
    [t,yHeun2] = explizitRK(F, t0, y0, T, n, 'Heun2');
    [t,yKutta3] = explizitRK(F, t0, y0, T, n, 'Kutta3');
    [t,yRK4] = explizitRK(F, t0, y0, T, n, 'RK4');
    
    ind = find(t == 3); % Approximation fuer t =3 (ind = Spaltenindex)
    % Alternativ in diesem Fall: ind = length(t);
    
    %2-Norm des Fehlers zwischen Approximation & analytische Lsg
    err(i-1,1) = norm(y_exact-yEuler(:,ind));
    err(i-1,2) = norm(y_exact-yHeun2(:,ind));
    err(i-1,3) = norm(y_exact-yKutta3(:,ind));
    err(i-1,4) = norm(y_exact-yRK4(:,ind));
end

% loglog-Plots (beide Axen logarithmisch skaliert)
loglog(h,err(:,1),h,err(:,2),h,err(:,3),h,err(:,4),h,h,'k--',h,h.^2,'k--',h,h.^3,'k--',h,h.^4,'k--')
title('loglog-Fehlerplot')
legend('Euler','Heun2','Kutta3','RK4','h, h^2, h^3, h^4')
xlabel('Schrittweite h')
ylabel('Fehler')