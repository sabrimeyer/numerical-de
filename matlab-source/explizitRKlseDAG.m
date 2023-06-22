function [t, y, z] = explizitRKlseDAG(A, B, C, D, v, w, t0, y0, T, n, methode)
    factor = B*D^-1; % Matrixprodukt
    dy = @(t,y) A*y-factor*(C*y+w(t))+v(t); % Differentialgleichung als Function-Handle

    [t,y] = explizitRK(dy, t0, y0, T, n, methode); % Approximation der Lsg
    m = length(y);
    
    z = zeros(length(D),m); % Spaltenweise Berechnung von (*) der Vektoren z_i
    for i =1:m
        z(:,i) = -inv(D)*(C*y(:,i) + w(t(i))); %(*) z = -D^-1*(C*y+w(t));
    end
end
