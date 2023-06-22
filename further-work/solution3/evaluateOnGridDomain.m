function fv = evaluateOnGridDomain(fh,G)
    nd = length(G);
    fv = zeros(nd,1);
    
    for k = 1:nd
        fv(k) = fh(G(1,k),G(2,k));
    end
    
end

