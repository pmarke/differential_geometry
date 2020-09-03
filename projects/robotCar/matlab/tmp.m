X1 = rand(6,1);
Y1 = rand(6,1);

x1 = se3_wedge(X1);
y1 = se3_wedge(Y1);

XX = expm(x1);
YY = expm(y1);

