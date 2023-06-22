function gv = evaluateOnGridBoundary(gh, B)
    nb = length(B);
    gv = zeros(nb,1);
    
    for k = 1:nb
        gv(k) = gh(B(1,k),B(2,k));
    end
    
end

