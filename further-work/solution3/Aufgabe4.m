%% Testbeispiel
% ------------------------------------------------------------------------
clc
clear
close all

u = @(x,y) abs(1./2.*x+y).*(1./2.*x+y) + 4.*x.^2.*y;
fh = @(x,y) -5/2*(.5*x+y>0) + 5/2*(.5*x+y<0) - 8*y;
% u(x,y) = g(x,y) auf dem Rand des Gebietes
gh = @(x,y) abs(1/2*x+y)*(1/2*x+y) + 4*x.^2*y;

n = 202;
p = 11;

[phi, c, x1, y1] = domains('face');
xs = x1(1):1/n:x1(2);
ys = y1(1):1/n:y1(2);

% xs = [-4:1/n:4];
% ys = [-4:1/n:4];

% phi = @(x,y) -x.^2 + -y.^2 + 0.25;
% c = 0;

[C, H, G, B] = makeGrid(phi, c, xs, ys, p);
[A, Ab] = shortleyWeller(C, H);
fv = evaluateOnGridDomain(fh,G);
gv = evaluateOnGridBoundary(gh, B); % u(x) = g(x) on boundary

uh = A\(fv - Ab*gv); % backslash-solving

figure

% Solution

subplot(2,2,1)
title(['Solution, n = ', num2str(n), ', p = ', num2str(p)])
xlabel('x')
ylabel('y')
zlabel('u')
[F, P] = makeMesh(C, H, G, B);
plotMeshFunction(F, P, uh, gv);

% Convergence
it = 11; %9
err = zeros(1,it);

for n = 3:it
    xs = x1(1):1/(2^n):x1(2);
    ys = y1(1):1/(2^n):y1(2);
    tic
    [C, H, G, B] = makeGrid(phi, c, xs, ys, p);
    Tgrid(n) = toc;
    tic
    [A, Ab] = shortleyWeller(C, H);
    Tshortley(n) = toc;
    tic
    fv = evaluateOnGridDomain(fh,G);
    gv = evaluateOnGridBoundary(gh, B); % u(x) = g(x) on boundary
    Teval(n) = toc;
    tic
    uh = A\(fv - Ab*gv); % backslash-solver
    Tsolver(n) = toc;
    err(n) = max(abs(uh' - u(G(1,:),G(2,:)))); % Fehler
    x(n) = length(G); % Freiheitsgrade
    h(n) = max(H, [], 'all');
end

subplot(2,2,2)
loglog(x, err , '-s', x, h.^2, '--', x, h, '--k')
legend('err', 'h^2', 'h')
title('Convergence')
xlabel('degrees of freedom, n')
ylabel('error, max-norm')

subplot(2,2,3)
loglog(x, Tgrid , '-s', x, Tshortley, '-s', x, Teval, '-s', x, Tsolver, '-s', x, 1e-5*x, '--k')
legend('makeGrid', 'shortleyWeller', 'evaluate', 'solver', 'n', 'Location', 'northwest')
title('Timings')
xlabel('degrees of freedom, n')
ylabel('time, s')

subplot(2,2,4)
loglog(Tshortley, err, '-s')
title('Work versus Error')
xlabel('time, s')
ylabel('error, max-norm')