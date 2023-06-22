function [t, y, z] = tiefpassRC(R, Cap, Ui, t0, y0, T, n, methode)
    A = 0;
    B = [0, Cap^-1];
    v = @(t) 0;
    C = [0; 1];
    D = [1 -R; 1 0];
    w = @(t) [0; -Ui(t)];
    
    [t, y, z] = explizitRKlseDAG(A, B, C, D, v, w, t0, y0, T, n, methode);
end

