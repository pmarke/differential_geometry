function gammaH= christoffelGeneralized(X,gammaG,G)

[~,c] = size(X);

gammaH = sym(zeros(c,c,c));

for delta=1:c
    n = (1/innerG(X(:,delta),X(:,delta),G));
    for alpha = 1:c
        for beta = 1:c
            A = affineConnections(X(:,alpha),X(:,beta),gammaG);
            gammaH(delta,alpha,beta) = n*innerG(A,X(:,delta),G);
        end
    end
    
end

gammaH= simplify(gammaH);


end

