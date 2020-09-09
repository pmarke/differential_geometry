


x = 0.2;





[z,fval] = fmincon(@fun,x,1,1)





%% Derivative test

rng('default')

P = eye(4)*2;
xh = [1;1;3;2];
x = [2;3;0;0];
C = [1 0 0 0; 0 1 0 0];
R = eye(2)*0.1;
r = chol(R);

v1 = r*randn(2,1);
v2 = r*randn(2,1);

y1 = C*x +v1;
y2 = C*x +v2;

z1 = y1 - C*xh;
z2 = y2 - C*xh;

Z = C'*inv(R)*z1 + C'*inv(R)*z2;

Pi = inv(P) + C'*inv(R)*C*2
xn = xh + inv(Pi)*Z


inv(Pi)*inv(P)*xh + inv(Pi)*C'*inv(R)*(y1+y2)


