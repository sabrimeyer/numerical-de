%%Test-Probleme für AUFGABE 3
% ----------------------------------------
%% Exponential-Funktion
clc
clear
close all

ya = 1;
f = @(x,y) y;
n = 15;
m = 20;

% [x,y] = odeABM2(f, 0, 10, ya, n, m);
[x,y] = odeABM3(f, 0, 10, ya, n, m);
plot(x,y,x,exp(x))
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
odeABM2(f, a, b, ya, n, m)
plot(x,y)