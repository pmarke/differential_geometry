function G = calculateMetricTensor(KE)
% Calculates the metric tensor from the kinetic energy.

syms x y z th phi Jw Jb mw mb l1 l2 'real'
q = [x y th phi];

G = sym(zeros(4,4));

for ii=1:4
    for jj = 1:4
        G(ii,jj) = diff(diff(KE,q(ii)),q(jj));
    end    
end
end

