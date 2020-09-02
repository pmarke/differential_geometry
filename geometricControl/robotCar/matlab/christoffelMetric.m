function gamma = christoffelMetric(G)
% Computes the christoffel symbols for the metric affine connection associated with G

syms x y z th phi Jw Jb mw mb l1 l2 mt 'real'

q = [x y th phi];

d = length(G);

gamma = sym(zeros(d,d,d));

for kk = 1:d
    for ii = 1:d
        for jj = 1:d
            tmp = 0;
            for ll = 1:d
                tmp = tmp + G(kk,ll)*(diff(G(ii,ll),q(jj))+diff(G(jj,ll),q(ii))-diff(G(ii,jj),q(ll)));
            end
            gamma(kk,ii,jj) =  (1/2)*tmp;
        end
    end
end

gamma = simplify(gamma);

end

