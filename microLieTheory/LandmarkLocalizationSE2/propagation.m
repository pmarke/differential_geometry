x_t = [1; 2; 0.1];
xh_t = [0.1; 0.2; 0.15];

u = [0.1; -0.2; 0.2];
dt = 0.01;
Q = diag([0.01,0.01,0.01]);
w = sqrt(Q)*randn(3,1);

x = Exp(x_t);
xh = Exp(xh_t);

e = Log(inv(xh)*x);

x_k = x*Exp((u+w*0)*dt);
xh_k = xh*Exp(u*dt);

e_ka = Log(inv(xh_k)*x_k)