%%Test-Probleme für AUFGABE 1
% ----------------------------------------
%% Exponential-Funktion
clc
clear
close all

ya = 1;
f = @(x, y) y;
% fy = 1;
n = 10;

[x, y] = odeAB2(f, 0, 10, ya, n);
[x1,y1] = odeAB3(f, 0, 10, ya, n);

alpha2 = [0 1/2 1/2 1];
beta2 = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
gamma2 = [1/6 1/3 1/3 1/6];
[x2, y2] = explizitRK(f, 0, 10, ya, n, alpha2, beta2, gamma2);

plot(x,y,x,exp(x),x1,y1,x2,y2)
legend('AB2','analytic','AB3', 'RK')

%% Lokta-Volterra-Modell
clc
clear
close all

alpha = 100;
f = @(x,y) [alpha*y(1)*(1 - y(2)), y(2)*(y(1) - 1)]';
fy = @(x,y) [alpha, -alpha*y(1); y(2), y(1)];
a = 0;
b = 5;
ya = [3;1];
sigma = 0.5;
tol = 1e-12;
maxiter = 42;
n=2^5;

%[x, y] = odeAB3(f, 0, 5, ya, n);
[x,y] = odeAM3(f, fy, a, b, ya, n, sigma, tol, maxiter);
plot(x,y)

%% Van-der-Pol Oszillator
clc
clear
close all

ya = [0, 2]';
n = 2^13;
c = 1/0.03;
f = @(x,y) [c*(y(2) - y(1).^3/3 + y(1)), -y(1)]';
% fy = @(x,y) [c*(1 - y(1).^2), c; -1, 0];
% [x, y] = odeAB2(f, 0, 10, ya, n);
[x, y] = odeAB3(f, 0, 10, ya, n);

alpha2 = [0 1/2 1/2 1];
beta2 = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
gamma2 = [1/6 1/3 1/3 1/6];
[x2, y2] = explizitRK(f, 0, 10, ya, n, alpha2, beta2, gamma2);

plot(x,y)