%% AUFGABE 2
%---------------------------------------------------------
%Beispiel y' = y
clc
clear all
close all

ya = 1;
f = @(x,y) y;
fy = @(x,y) 1;
n = 2^6;
theta = 0.5;
sigma = 0.5;
tol = 1e-12;
maxiter = 42;

[x,y] = thetaSchema(f, fy, 1, 3, ya, n, theta, sigma, tol, maxiter);

fehler = norm(exp(x-1) - y)
% plot(x,y,x,exp(x-1))