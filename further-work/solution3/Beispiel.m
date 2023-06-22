%% Zu Aufgabe 1
clc
clear
close all

n = 1;
R = 10;
xs = [-R:1/n:R];
ys = [-R:1/n:R];

phi = @(x,y) -x.^2 + - y.^2 + 5;
c = 2;

[C, H, G, B] = makeGrid(phi, c, xs, ys, 10);

figure();
plotGrid(C, H, G, B);

figure();
plotGridMesh(C, H, G, B);

%% Zu Aufgabe 2

clc
clear
close all

G = [.5 , .2 , .5 , .7 , .2 , .5 , .7 , 1;
     .1 , .3 , .3 , .3 , .5 , .5 , .5 , .5];
B = [.28 , .63 , .5 , .08 , .2 , .92 , .7 , .13 , .2 , .5 , .7 , 1.05 , 1.0 , 1.0;
     .1 , .1 , .05 , .3 , .15 , .3 , .14 , .5 , .57 , .67 , .63 , .5 , .38 , .58];
C = [-1 , -4 ,  2 ,  3 , -8 ,   5 ,   6 ,   7;
     -2 ,  3 ,  4 , -6 ,  6 ,   7 ,   8 , -12;
     -3 , -5 ,  1 , -7 ,  2 ,   3 ,   4 , -13;
      3 ,  5 ,  6 ,  7 , -9 , -10 , -11 , -14];
H = [.22 , .12 , .3 , .2  , .07 , .3  , .2  , .3;
     .13 , .3  , .2 , .22 , .3  , .2  , .3  , .05;
     .05 , .15 , .2 , .16 , .2  , .2  , .2  , .12;
     .2  , .2  , .2 , .2  , .07 , .17 , .13 , .08];
 
[A, Ab] = shortleyWeller(C, H);

%% Zu Aufgabe 3

clc
clear
close all

n = 25;
R = 5;
xs = [-R:1/n:R];
ys = [-R:1/n:R];

phi = @(x,y) -x.^2 + -y.^2 + 5;
c = 0;

[C, H, G, B] = makeGrid(phi, c, xs, ys, 10);
 
[A, Ab] = shortleyWeller(C, H);

fh = @(x,y) x.^2 + y.^2;
gh = @(x,y) 0;

fv = evaluateOnGridDomain(fh,G);
gv = evaluateOnGridBoundary(gh, B); % u(x) = g(x) on boundary

u = A\(fv - Ab*gv); % backslash-solving

% [F, P] = makeMesh(C, H, G, B);
% plotMeshFunction(F, P, u, gv);
