%% AUFGABE 3
%---------------------------------------------------------
%Beispiel y' = y
clc
clear all
close all

ya = 1;
f = @(x, y) y;
n = 100;

% alpha = 0;
% beta = 0;
% gamma = 1;

% alpha = [0 1];
% beta = [0 0; 1 0];
% gamma = [0.5 0.5];

alpha = [0 0.5];
beta = [0 0; 0.5 0];
gamma = [0 1];


[x, y] = explizitRK(f, 1, 3, ya, n, alpha, beta, gamma);

fehler = norm(exp(x-1) - y)
% plot(x,y,x,exp(x-1))