%%Test-Probleme für AUFGABE 2
% ----------------------------------------
%% Exponential-Funktion
clc
clear
close all

ya = 1;
f = @(x,y) y;
fy = @(x,y) 1;
n = 20;
theta = 0.5;
sigma = 0.5;
tol = 1e-12;
maxiter = 42;

[x, y] = odeAM2(f, fy, 0, 10, ya, n, sigma, tol, maxiter);
% [x, y] = odeAM3(f, fy, 0, 10, ya, n, sigma, tol, maxiter);
% [x,y] = odeBDF2(f, fy, 0, 10, ya, n, sigma, tol, maxiter);
[xx,yy] = odeABM2(f, 0, 10, ya, n, 20);
% [x,y] = odeBDF3(f, fy, 0, 10, ya, n, sigma, tol, maxiter);
plot(x,y,x,exp(x),xx, yy)
legend('numerical', 'analytic')

%% Van-der-Pol Oszillator
clc
clear
close all

theta = 0.5;
sigma = 0.5;
tol = 1e-12;
maxiter = 42;

ya = [0, 2]';
c = 1/0.03;
f = @(x,y) [c*(y(2) - y(1).^3/3 + y(1)), -y(1)]';
fy = @(x,y) [c*(1 - y(1).^2), c; -1, 0];
n = 2^9;
% [x,y] = odeAM2(f, fy, 0, 10, ya, n, 0.5, 1e-6, 1e4);
% [x, y] = odeAM3(f, fy, 0, 10, ya, n, sigma, tol, maxiter);
% [x,y] = odeBDF2(f, fy, 0, 10, ya, n, sigma, tol, maxiter);
[x,y] = odeBDF3(f, fy, 0, 10, ya, n, sigma, tol, maxiter);
plot(x,y)