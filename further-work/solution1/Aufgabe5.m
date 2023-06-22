%% AUFGABE 5: (Beispiel I: Testproblem)
% -------------------------------------------------------------------------
%% (a) thetaSchema
clear all
close all
clc

sigma = 0.5;
tol = 1e-12;
maxiter = 42;
theta = [0 0.49 0.499 0.4999 0.5 0.5001 0.501 0.51 1];

errorMatrix = dataThetaSchema(1, theta, sigma, tol, maxiter);
errorMatrix2 = dataThetaSchema(-1, theta, sigma, tol, maxiter);

n = 2.^(3:13);
figure
sgtitle('Fehler gegen Iterationen n (für Lambda = 1)')
for j = 1:9
    subplot(3,3,j)
    loglog(n, errorMatrix(:,j), n, n.^(-1), 'k-.', n, n.^(-2), 'k-.')
    title(num2str(theta(j)))
end

figure
sgtitle('Fehler gegen Iterationen n (für Lambda = -1)')
for j = 1:9
    subplot(3,3,j)
    loglog(n, errorMatrix2(:,j), n, n.^(-1), 'k-.', n, n.^(-2), 'k-.')
    title(num2str(theta(j)))
end


%% (b) explizitRK
clc
clear all
close all

errorMatrix = dataExplizitRK(1, 'Euler');
errorMatrix2 = dataExplizitRK(1, 'Heun');
errorMatrix3 = dataExplizitRK(1, 'Klassisch');

errorMatrix4 = dataExplizitRK(-1, 'Euler');
errorMatrix5 = dataExplizitRK(-1, 'Heun');
errorMatrix6 = dataExplizitRK(-1, 'Klassisch');

n = 2.^(3:13);
figure
sgtitle('Fehler gegen Iterationen n')
subplot(1,2,1)
loglog(n, errorMatrix, n, errorMatrix2, n, errorMatrix3, n, n.^(-1), 'k-.', n, n.^(-2), 'k-.', n, n.^(-4), 'k-.')
title('Lambda = 1')

subplot(1,2,2)
loglog(n, errorMatrix4, n, errorMatrix5, n, errorMatrix6, n, n.^(-1), 'k-.', n, n.^(-2), 'k-.', n, n.^(-4), 'k-.')
title('Lambda = -1')
legend('Euler', 'Heun', 'Klassisch', 'Raten', 'location', 'southwest')
legend('boxoff')

%% (c) explizitRKadapt
clc
clear all
close all

m = 12; % zu langsam für m = 20 !

errorMatrix = dataExplizitRKadapt(1, 'RKF2(3)', m);
errorMatrix2 = dataExplizitRKadapt(1, 'DOPRI', m);

errorMatrix3 = dataExplizitRKadapt(-1, 'RKF2(3)', m);
errorMatrix4 = dataExplizitRKadapt(-1, 'DOPRI', m);

axis = 2.^(-1:-1:-m);
figure
sgtitle('Fehler gegen Genauigkeit epsilon')
subplot(1,2,1)
loglog(axis, errorMatrix, axis, errorMatrix2, axis, axis, 'r-.', axis, axis.^(2/3), 'k-.', axis, axis.^(4/5), 'k-.')
title('Lambda = 1')

subplot(1,2,2)
loglog(axis, errorMatrix3, axis, errorMatrix4, axis, axis, 'r-.', axis, axis.^(2/3), 'k-.', axis, axis.^(4/5), 'k-.')
title('Lambda = -1')
legend('RKF2(3)', 'DOPRI', 'linear', 'location', 'southeast')
legend('boxoff')

%% Funktionen
function errorMatrix = dataThetaSchema(lambda, theta, sigma, tol, maxiter)
errorMatrix = zeros(11, 9);
f = @(x,y) lambda*y;
fy = @(x,y) lambda;

    for i = 1:11
        n = 2^(i + 2);
        for j = 1:length(theta)
            [x, y] = thetaSchema(f, fy, 1, 3, 1, n, theta(j), sigma, tol, maxiter);
            errorMatrix(i, j) = abs(y(length(x)) - exp(lambda*(3 - 1)));
        end
    end
end

function errorMatrix = dataExplizitRK(lambda, verfahren)
errorMatrix = zeros(11,1);
f = @(x,y) lambda*y;

switch verfahren
    case 'Euler'
        alpha = 0;
        beta = 0;
        gamma = 1;
        
    case 'Heun'
        alpha = [0 1];
        gamma = [1/2 1/2];
        beta = [0 0; 1 0];
        
    case 'Klassisch'
        alpha = [0 1/2 1/2 1];
        beta = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
        gamma = [1/6 1/3 1/3 1/6];
end

    for i = 1:11
        n = 2^(i + 3);
        [x, y] = explizitRK(f, 1, 3, 1, n, alpha, beta, gamma);
        m = length(x);
        errorMatrix(i) = abs(y(m) - exp(lambda*(3 - 1)));
    end
end

function errorMatrix = dataExplizitRKadapt(lambda, verfahren, m)
errorMatrix = zeros(m,1);
f = @(x,y) lambda*y;

switch verfahren
    case 'RKF2(3)'
        alpha = [0 1 1/2];
        beta = [0 0 0; 1 0 0; 1/4 1/4 0];
        gamma = [1/2 1/2 0];
        gammah = [1/6 1/6 4/6];
        p = 2;
        
    case 'DOPRI'
        alpha = [0 1/5 3/10 4/5 8/9 1 1];
        beta = [0 0 0 0 0 0 0; 1/5 0 0 0 0 0 0; 3/40 9/40 0 0 0 0 0; 44/45 -56/15 32/9 0 0 0 0; 19372/6561 -25360/2187 64448/6561 -212/729 0 0 0; 9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0; 35/384 0 500/1113 125/192 -2187/6784 11/84 0];
        gamma = [35/384 0 500/1113 125/192 -2187/6784 11/84 0];
        gammah = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40];
        p = 4;
        
end

    for i = 1:m
        epsilon = 2^(-i);
        [x, y] = explizitRKadapt(f, 1, 3, 1, 0.5, alpha, beta, gamma, gammah, p, 0.9, epsilon);
        m = length(x);
        errorMatrix(i) = abs(y(m) - exp(lambda*(x(end) - 1)));
    end
end
