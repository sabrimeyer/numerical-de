function [A, Ab] = shortleyWeller(C, H)
    nd = length(C);
    ai = [1 1 1 1 1]';
    aj = [1, C(1,1), C(2,1), C(3,1), C(4,1)]';
    av = [2/(H(1,1)*H(2,1))+2/(H(3,1)*H(4,1)); -2/(H(1,1)*(H(1,1)+H(2,1)));
        -2/(H(2,1)*(H(1,1)+H(2,1))); -2/(H(3,1)*(H(3,1)+H(4,1))); -2/(H(4,1)*(H(3,1)+H(4,1)))];
    
    % Spaltenvektoren initialisieren
    hi = @(k) [k k k k k]';
    hj = @(k) [k C(1,k) C(2,k) C(3,k) C(4,k)]';
    hv = @(k) [2/(H(1,k)*H(2,k))+2/(H(3,k)*H(4,k)); -2/(H(1,k)*(H(1,k)+H(2,k)));
        -2/(H(2,k)*(H(1,k)+H(2,k))); -2/(H(3,k)*(H(3,k)+H(4,k))); -2/(H(4,k)*(H(3,k)+H(4,k)))];;
    
    for s = 2:nd
        ai = [ai; hi(s)];
        aj = [aj; hj(s)];
        av = [av; hv(s)];
    end
    
    % Spaltenvektoren modifizieren
    ai_hat = ai;
    aj_hat = aj;
    av_hat = av;
    ai_til = ai;
    aj_til = aj;
    av_til = av;
    
    ai_hat(aj<0) = [];
    aj_hat(aj<0) = [];
    av_hat(aj<0) = [];
    ai_til(aj>0) = [];
    aj_til(aj>0) = [];
    av_til(aj>0) = [];
    
    % Matrzen aufbauen
    A = sparse(ai_hat, aj_hat, av_hat);
    Ab = sparse(ai_til, -aj_til, av_til);
    
end
