
%% AUFGABE 1
%-------------------------------------------------------
clc
clear all
close all

g = @(x) x^2 - 2;
gp = @(x) 2*x;
z0 = 1;
newtonIterationSWS(g, gp, z0, 0.5, 0.0001, 1000)

%% Nichtlineares Gleichungssystem
clc
clear all
close all

g = @(x) [x(1).^2 + x(2).^2 - 9; x(1) + x(2) - 3];
gp = @(x) [2*x(1), 2*x(2); 1, 1];
z0 = [100, 1]';
z = newtonIterationSWS(g, gp, z0, 0.5, 0.0001, 1000)
