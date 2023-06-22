%% AUFGABE 8: (Beispiel IV: SEIRS-Modell)
% -------------------------------------------------------------------------
% Fall 0
clc
clear all
close all

% rEE = 4;
% rEI = 6;

[x, y, x2, y2, x3, y3] = getDataFall0(0, 12);

figure

subplot(3,3,1)
area(x', y')
title('thetaSchema')
xlim([0,12])
legend('R', 'S', 'E', 'I', 'location', 'southeast')
ylabel('Fall 0')

subplot(3,3,2)
area(x2', y2')
xlim([0,12])
title('explizitRK (Klassisch)')

subplot(3,3,3)
area(x3', y3')
xlim([0,12])
title('explizitRKadapt (DP)')

[x, y, x2, y2, x3, y3] = getDataFall1(0, 12);

subplot(3,3,4)
area(x', y')
xlim([0,12])
ylabel('Fall 1')

subplot(3,3,5)
area(x2', y2')
xlim([0,12])

subplot(3,3,6)
area(x3', y3')
xlim([0,12])

[x, y, x2, y2, x3, y3] = getDataFall2(0, 12);

subplot(3,3,7)
area(x', y')
xlim([0,12])
ylabel('Fall 2')

subplot(3,3,8)
area(x2', y2')
xlim([0,12])

subplot(3,3,9)
area(x3', y3')
xlim([0,12])


%% Funktionen
function [x, y, x2, y2, x3, y3] = getDataFall0(a, b)

    rR = 2;
    rS = 0.3;
    rI = 3;

    f = @(x,y) [rR*y(4) - rS*y(1), rS*y(1) - rEE(0)*y(2).*y(3) - rEI(0)*y(2).*y(4), rEE(0)*y(2).*y(3) + rEI(0)*y(2).*y(4) - rI*y(3), rI*y(3) - rR*y(4)]';
    fy = @(x,y) [-rS, 0, 0, rR; rS, -rEE(0).*y(3) - rEI(0).*y(4), -rEE(0).*y(2), -rEI(0).*y(2); 0, rEE(0).*y(3) + rEI(0).*y(4), rEE(0).*y(2) - rI, rEI(0).*y(2); 0, 0, rI, -rR];
    ya = [0, 0.999999, 0.000001, 0]';
    [x, y] = thetaSchema(f, fy, a, b, ya, 2^8, 0.6, 0.5, 1e-12, 42);

    alpha2 = [0 1/2 1/2 1];
    beta2 = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
    gamma2 = [1/6 1/3 1/3 1/6];
    [x2, y2] = explizitRK(f, a, b, ya, 2^8, alpha2, beta2, gamma2);

    alpha3 = [0 1/5 3/10 4/5 8/9 1 1];
    beta3 = [0 0 0 0 0 0 0; 1/5 0 0 0 0 0 0; 3/40 9/40 0 0 0 0 0; 44/45 -56/15 32/9 0 0 0 0; 19372/6561 -25360/2187 64448/6561 -212/729 0 0 0; 9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0; 35/384 0 500/1113 125/192 -2187/6784 11/84 0];
    gamma3 = [35/384 0 500/1113 125/192 -2187/6784 11/84 0];
    gammah = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40];
    p = 4;
    epsilon = 1e-5; % 1e-5 zu langsam!
    [x3, y3] = explizitRKadapt(f, a, b, ya, 0.1, alpha3, beta3, gamma3, gammah, p, 0.9, epsilon);
    
end

function [x, y, x2, y2, x3, y3] = getDataFall1(a, b)

    rR = 2;
    rS = 0.3;
    rI = 3;

    f = @(x,y) [rR*y(4) - rS*y(1), rS*y(1) - rEE(0)*y(2).*y(3) - rEI(x)*y(2).*y(4), rEE(0)*y(2).*y(3) + rEI(x)*y(2).*y(4) - rI*y(3), rI*y(3) - rR*y(4)]';
    fy = @(x,y) [-rS, 0, 0, rR; rS, -rEE(0).*y(3) - rEI(x).*y(4), -rEE(0).*y(2), -rEI(x).*y(2); 0, rEE(0).*y(3) + rEI(x).*y(4), rEE(0).*y(2) - rI, rEI(x).*y(2); 0, 0, rI, -rR];
    ya = [0, 0.999999, 0.000001, 0]';
    [x, y] = thetaSchema(f, fy, a, b, ya, 2^8, 0.6, 0.5, 1e-12, 42);

    alpha2 = [0 1/2 1/2 1];
    beta2 = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
    gamma2 = [1/6 1/3 1/3 1/6];
    [x2, y2] = explizitRK(f, a, b, ya, 2^8, alpha2, beta2, gamma2);

    alpha3 = [0 1/5 3/10 4/5 8/9 1 1];
    beta3 = [0 0 0 0 0 0 0; 1/5 0 0 0 0 0 0; 3/40 9/40 0 0 0 0 0; 44/45 -56/15 32/9 0 0 0 0; 19372/6561 -25360/2187 64448/6561 -212/729 0 0 0; 9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0; 35/384 0 500/1113 125/192 -2187/6784 11/84 0];
    gamma3 = [35/384 0 500/1113 125/192 -2187/6784 11/84 0];
    gammah = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40];
    p = 4;
    epsilon = 1e-5; % 1e-5 zu langsam!
    [x3, y3] = explizitRKadapt(f, a, b, ya, 0.1, alpha3, beta3, gamma3, gammah, p, 0.9, epsilon);
    
end

function [x, y, x2, y2, x3, y3] = getDataFall2(a, b)

    rR = 2;
    rS = 0.3;
    rI = 3;

    f = @(x,y) [rR*y(4) - rS*y(1), rS*y(1) - rEE(x)*y(2).*y(3) - rEI(x)*y(2).*y(4), rEE(x)*y(2).*y(3) + rEI(x)*y(2).*y(4) - rI*y(3), rI*y(3) - rR*y(4)]';
    fy = @(x,y) [-rS, 0, 0, rR; rS, -rEE(x).*y(3) - rEI(x).*y(4), -rEE(x).*y(2), -rEI(x).*y(2); 0, rEE(x).*y(3) + rEI(x).*y(4), rEE(x).*y(2) - rI, rEI(x).*y(2); 0, 0, rI, -rR];
    ya = [0, 0.999999, 0.000001, 0]';
    [x, y] = thetaSchema(f, fy, a, b, ya, 2^8, 0.6, 0.5, 1e-12, 42);

    alpha2 = [0 1/2 1/2 1];
    beta2 = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
    gamma2 = [1/6 1/3 1/3 1/6];
    [x2, y2] = explizitRK(f, a, b, ya, 2^8, alpha2, beta2, gamma2);

    alpha3 = [0 1/5 3/10 4/5 8/9 1 1];
    beta3 = [0 0 0 0 0 0 0; 1/5 0 0 0 0 0 0; 3/40 9/40 0 0 0 0 0; 44/45 -56/15 32/9 0 0 0 0; 19372/6561 -25360/2187 64448/6561 -212/729 0 0 0; 9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0; 35/384 0 500/1113 125/192 -2187/6784 11/84 0];
    gamma3 = [35/384 0 500/1113 125/192 -2187/6784 11/84 0];
    gammah = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40];
    p = 4;
    epsilon = 1e-5; % 1e-5 zu langsam!
    [x3, y3] = explizitRKadapt(f, a, b, ya, 0.1, alpha3, beta3, gamma3, gammah, p, 0.9, epsilon);
    
end

function val = rEE(x)
    if x < 3
        val = 4;
    else
        val = 2;
    end
end

function val = rEI(x)
    if x < 2
        val = 6;
    else
        val = 0.5;
    end
end
