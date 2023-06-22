%% Testbeispiel
% ------------------------------------------------------------------------
clc
clear
close all

n = 50;
gh = @(x,y) 0;
%fh = @(x,y) 0;

[phi, c, xl, yl] = domains('ellipsering');
xs = xl(1):1/n:xl(2);
ys = yl(1):1/n:yl(2);

[C, H, G, B] = makeGrid(phi, c, xs, ys, 10);

[A, Ab] = shortleyWeller(C, H);
[V, D] = eigs(A,10,'smallestabs');
%A = sparse(A-100*eye(size(A)));
lambda = diag(D);

figure

for it = 1:10
    %A = sparse(A-lambda(it)*eye(size(A)));
    %fv = evaluateOnGridDomain(fh,G);
    fv = V(:,it);
    gv = evaluateOnGridBoundary(gh, B); % u(x) = g(x) on boundary

    u = A\fv; % backslash-solving
    u = 0.2*u/max(u); % normalize data

    subplot(2,5,it)
    title(num2str(lambda(it)))
    [F, P] = makeMesh(C, H, G, B);
    xlim([-.5,.5])
    ylim([-.5,.5])
    zlim([-.25,.25])
    plotMeshFunction(F, P, u, gv);
end