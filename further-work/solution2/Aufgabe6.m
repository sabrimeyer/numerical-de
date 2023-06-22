%% AUFGABE 4 (Beispiel III: Van-der-Pol-Oszillator)
% -------------------------------------------------------------------------
clc
clear
close all

c = 1/0.03;
f = @(x,y) [c*(y(2) - y(1).^3/3 + y(1)), -y(1)]';
fy = @(x,y) [c*(1 - y(1).^2), c; -1, 0];
a = 0;
b = 10;
ya = [0;2];
sigma = 0.5;
tol = 1e-12;
maxiter = 42;
m = 2;
exact = [1.556042079300293; -0.332027354872852];

for i = 3:13
    n = 2^i;
    [x,y] = odeAB2(f, a, b, ya, n);
    [x1,y1] = odeAB3(f, a, b, ya, n);
    [x2,y2] = odeAM2(f, fy, a, b, ya, n, sigma, tol, maxiter);
    [x3,y3] = odeAM3(f, fy, a, b, ya, n, sigma, tol, maxiter);
    [x4,y4] = odeBDF2(f, fy, a, b, ya, n, sigma, tol, maxiter);
    [x5,y5] = odeBDF3(f, fy, a, b, ya, n, sigma, tol, maxiter);
    [x6,y6] = odeABM2(f, a, b, ya, n, m);
    [x7,y7] = odeABM3(f, a, b, ya, n, m);
    errorAB2(i-2) = norm(y(:,end) - exact);
    errorAB3(i-2) = norm(y1(:,end) - exact);
    errorAM2(i-2) = norm(y2(:,end) - exact);
    errorAM3(i-2) = norm(y3(:,end) - exact);
    errorBDF2(i-2) = norm(y4(:,end) - exact);
    errorBDF3(i-2) = norm(y5(:,end) - exact);
    errorABM2(i-2) = norm(y6(:,end) - exact);
    errorABM3(i-2) = norm(y7(:,end) - exact);
end

n = 2.^(3:13);
c = 1e5;

figure

sgtitle('Mehrschrittverfahren, Van-der-Pol: \epsilon = 0.03')

subplot(2,2,1)
loglog(n, errorAB2, '-s', n, errorAB3, '-s', n, c*n.^-2, '-.', n, c*n.^-3, '-.', n, c*n.^-4, '-.')
legend('AB2', 'AB3', 'n^{-2}', 'n^{-3}', 'n^{-4}', 'Location', 'southwest')
title('Adams-Bashforth')
xlabel('Anzahl Zeitschritte n')
ylabel('Fehler bei x=3')

subplot(2,2,2)
loglog(n, errorAM2, '-s', n, errorAM3, '-s', n, c*n.^-2, '-.', n, c*n.^-3, '-.', n, c*n.^-4, '-.')
legend('AM2', 'AM3', 'n^{-2}', 'n^{-3}', 'n^{-4}', 'Location', 'southwest')
title('Adams-Moulton')
xlabel('Anzahl Zeitschritte n')
ylabel('Fehler bei x=3')

subplot(2,2,3)
loglog(n, errorBDF2, '-s', n, errorBDF3, '-s', n, c*n.^-2, '-.', n, c*n.^-3, '-.', n, c*n.^-4, '-.')
legend('BDF2', 'BDF3', 'n^{-2}', 'n^{-3}', 'n^{-4}', 'Location', 'southwest')
title('BDF')
xlabel('Anzahl Zeitschritte n')
ylabel('Fehler bei x=3')

subplot(2,2,4)
loglog(n, errorABM2, '-s', n, errorABM3, '-s', n, c*n.^-2, '-.', n, c*n.^-3, '-.', n, c*n.^-4, '-.')
legend('ABM2', 'ABM3', 'n^{-2}', 'n^{-3}', 'n^{-4}', 'Location', 'southwest')
title('Adams-Bashforth-Moulton, m = 2')
xlabel('Anzahl Zeitschritte n')
ylabel('Fehler bei x=3')





