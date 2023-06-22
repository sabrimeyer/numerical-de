%% AUFGABE 6: Beispiel II (Lotka-Volterra-Modell)
% -------------------------------------------------------------------------
clc
clear all
close all

f = @(x,y) [10*y(1)*(1 - y(2)), y(2)*(y(1) - 1)]';
fy = @(x,y) [10 -10*y(1); y(2) y(1)]; %Jacobi-Matrix
ya = [3 1]';

[x, y] = thetaSchema(f, fy, 0, 5, ya, 2^6, 0.5001, 0.5, 1e-12, 42);

alpha2 = [0 1/2 1/2 1];
beta2 = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
gamma2 = [1/6 1/3 1/3 1/6];
[x2, y2] = explizitRK(f, 0, 5, ya, 2^6, alpha2, beta2, gamma2);

% alpha3 = [0 1/2 1/2 1 1];
% beta3 = [0 0 0 0 0; 1/2 0 0 0 0; 0 1/2 0 0 0; 0 0 1 0 0; 1/6 1/3 1/3 1/6 0];
% gamma3 = [1/6 1/3 1/3 0 1/6];
% gammah = [1/6 1/3 1/3 1/6 0];
% p = 3;

% alpha3 = [0 1 1/2];
% beta3 = [0 0 0; 1 0 0; 1/4 1/4 0];
% gamma3 = [1/2 1/2 0];
% gammah = [1/6 1/6 4/6];
% p = 2;

alpha3 = [0 1/5 3/10 4/5 8/9 1 1];
beta3 = [0 0 0 0 0 0 0; 1/5 0 0 0 0 0 0; 3/40 9/40 0 0 0 0 0; 44/45 -56/15 32/9 0 0 0 0; 19372/6561 -25360/2187 64448/6561 -212/729 0 0 0; 9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0; 35/384 0 500/1113 125/192 -2187/6784 11/84 0];
gamma3 = [35/384 0 500/1113 125/192 -2187/6784 11/84 0];
gammah = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40];
p = 4;
epsilon = 1e-2; % zu langsam für 1e-5 (und Punkte schwer erkennbar)
[x3, y3] = explizitRKadapt(f, 0, 5, ya, 1, alpha3, beta3, gamma3, gammah, p, 0.9, 1e-2);


figure

subplot(1,3,1)
plot(x, y(1,:), x, y(2,:), x, 0, 'r.')
title('thetaSchema')

subplot(1,3,2)
plot(x2, y2(1,:), x2, y2(2,:), x2, 0, 'r.')
title('explizitRK (Klassisch)')

subplot(1,3,3)
plot(x3, y3(1,:), x3, y3(2,:), x3, 0, 'r.')
title('explizitRKadapt (Dormand-Prince)')
legend('y(1)', 'y(2)', 'gridpoints', 'location', 'southeast')
xlim([0,5])
ylim([0,3])

