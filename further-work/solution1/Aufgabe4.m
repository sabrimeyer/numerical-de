%% AUFGABE 3
%---------------------------------------------------------
%Beispiel y' = y
clc
clear all
close all

ya = 1;
f = @(x, y) y;

% alpha = 0;
% beta = 0;
% gamma = 1;
% gammah = 0.1;

alpha = [0 1 0.5];
beta = [0 0 0; 1 0 0; 0.25 0.25 0];
gamma = [0.5 0.5 0];
gammah = [1/6 1/6 4/6];
p = 2;
h0 = 0.5;
epsilon = 0.001;

[x, y] = explizitRKadapt(f, 1, 3, ya, h0, alpha, beta, gamma, gammah, p, 0.9, epsilon);

fehler = norm(exp(x-1) - y)

% plot(x,y,x,exp(x-1))