function [C, H, G, B] = makeGrid(phi, c, xs, ys, p)
    r = length(xs);
    s = length(ys);
    
    
    % Berechne I
    I = zeros(r,s);
    nd = 1;
    for k = 1:s
        for j = 1:r
            if phi(xs(j),ys(k)) > c
                I(j,k) = nd;
                nd = nd + 1;
            end
        end
    end
   nd = nd - 1; % ESSENTIELLE KORREKTUR !!!
    
    % Berechne C,G
    C = zeros(4,nd);
    G = zeros(2,nd);
    for j = 1:r
        for k = 1:s
            if I(j,k) ~= 0
                C(1,I(j,k)) = I(j-1,k);
                C(2,I(j,k)) = I(j+1,k);
                C(3,I(j,k)) = I(j,k-1);
                C(4,I(j,k)) = I(j,k+1);
                
                G(1,I(j,k)) = xs(j);
                G(2,I(j,k)) = ys(k);
            end
        end
    end
    idx = C == 0;
    nb = sum(idx(:)); % zeros in C counted
    
    % Aktualisiere C & berechne B und H
    m = 0;
    B = zeros(2, nb);
    H = zeros(4,nd);
    for i = 1:nd
        for q = 1:4
            if C(q,i) == 0 % Randpunkt gefunden
                m = m + 1;
                C(q,i) = -m;
                if q < 3
                    B(1,m) = bisection(G(:,i),q, xs, ys, p, phi, c); % Bisektionsschritte für x
                    B(2,m) = G(2,i);
                    H(q,i) = norm(G(:,i) - B(:,m));
                else
                    B(1,m) = G(1,i);
                    B(2,m) = bisection(G(:,i),q, xs, ys, p, phi, c); % Bisektionsschritte für y
                    H(q,i) = norm(G(:,i) - B(:,m));
                end
            else
                H(q,i) = norm(G(:,i) - G(:,C(q,i)));
            end
        end 
    end
    
end

function b = bisection(pt, dir, xs, ys, p, phi, c) % point & direction
    switch dir
        
        case 1 % links
            out = xs(find(xs == pt(1)) - 1);
            b = pt(1);
            h = (out - pt(1))/p; % Bisektionsschrittweite
            while phi(b, pt(2)) > c && b > out
                b = b + h;
            end
            
        case 2 % rechts
            out = xs(find(xs == pt(1)) + 1);
            b = pt(1);
            h = (out - pt(1))/p; % Bisektionsschrittweite
            while phi(b, pt(2)) > c && b < out
                b = b + h;
            end
            
        case 3 % unten
            out = ys(find(ys == pt(2)) - 1);
            b = pt(2);
            h = (out - pt(2))/p; % Bisektionsschrittweite
            while phi(pt(1), b) > c && b > out
                b = b + h;
            end
            
        case 4 % oben
            out = ys(find(ys == pt(2)) + 1);
            b = pt(2);
            h = (out - pt(2))/p; % Bisektionsschrittweite
            while phi(pt(1), b) > c && b < out
                b = b + h;
            end
            
    end
end

