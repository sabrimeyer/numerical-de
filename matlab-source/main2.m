close all
clear all
clc

for k = 1:3 % n = 10000, 20000, 30000

%Initialisierung
n = k*10000; % Anzahl Datenpunkte: n+1
R = 0.5; % Widerstand
Cap = 20.23e-4; % Kapazitaet des Kondensators
t0 = 0; % linker Randpunkt
T = 2; % rechter Randpunkt
Ui = @(t) sin(2*pi*t.*(220*2.^t)); % zugefuehrte Spannung
y0 = Ui(0); % Anfangswert

[t, yEuler, z] = tiefpassRC(R, Cap, Ui, t0, y0, T, n, 'Euler');
[t, yHeun2, z] = tiefpassRC(R, Cap, Ui, t0, y0, T, n, 'Heun2');
[t, yKutta3, z] = tiefpassRC(R, Cap, Ui, t0, y0, T, n, 'Kutta3');
[t, yRK4, z] = tiefpassRC(R, Cap, Ui, t0, y0, T, n, 'RK4');


figure

subplot(2,2,1)
plot(t,Ui(t),t,yEuler)
title('Euler')
xlabel('Zeit t')
ylabel('Spannung U')
legend('U_i','U_2')


subplot(2,2,2)
plot(t,Ui(t),t,yHeun2)
title('Heun2')

subplot(2,2,3)
plot(t,Ui(t),t,yKutta3)
title('Kutta3')

subplot(2,2,4)
plot(t,Ui(t),t,yRK4)
title('RK4')

sgtitle(['n = ', num2str(k*10000)])
end


