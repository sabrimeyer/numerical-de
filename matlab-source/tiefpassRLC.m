function [t, y, z] = tiefpassRLC(R, L, Cap, Ui, t0, y0, T, n, methode)
    A = zeros(2);
    B = [0 0 0 Cap^-1; 0 L^-1 0 0];
    v = @(t) [0; 0];
    C = [0 0; 0 1; -1 0; 1 0];
    D = [1 0 -R 0; 0 0 -1 -1; 1 0 0 0; 0 1 0 0];
    w = @(t) [0; 0; 0; -Ui(t)];
    
    [t, y, z] = explizitRKlseDAG(A, B, C, D, v, w, t0, y0, T, n, methode);
end