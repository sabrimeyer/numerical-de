clc
clear all
close all

for k = 1:3 % n = 10000, 20000, 30000

n = k*10000;
R = 0.5;
L = 3.61e-4;
Cap = 3.61e-4;
t0 = 0;
T = 2;
Ui = @(t) sin(2*pi*t.*(220*2.^t));
y0 = [0 0];

[t, yEuler, z] = tiefpassRLC(R, L, Cap, Ui, t0, y0, T, n, 'Euler');
[t, yHeun2, z] = tiefpassRLC(R, L, Cap, Ui, t0, y0, T, n, 'Heun2');
[t, yKutta3, z] = tiefpassRLC(R, L, Cap, Ui, t0, y0, T, n, 'Kutta3');
[t, yRK4, z] = tiefpassRLC(R, L, Cap, Ui, t0, y0, T, n, 'RK4');

figure

subplot(2,2,1)
plot(t,Ui(t),t,yEuler(1,:)) %nur erste Zeile weil yEuler = [U3,I2]
title('Euler')
xlabel('Zeit t')
ylabel('Spannung U')
legend('U_i','U_3')


subplot(2,2,2)
plot(t,Ui(t),t,yHeun2(1,:))
title('Heun2')

subplot(2,2,3)
plot(t,Ui(t),t,yKutta3(1,:))
title('Kutta3')

subplot(2,2,4)
plot(t,Ui(t),t,yRK4(1,:))
title('RK4')

sgtitle(['n = ', num2str(k*10000)])

end
