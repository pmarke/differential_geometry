function M  = se3_adx(X)

p = X(1:3);
w = X(4:6);

M = [ssm(w) ssm(p); zeros(3,3) ssm(w)];

end