%
rng('default')
E0 = [eye(3) zeros(3,1);zeros(1,3) 1];
W0 = [SE3_dse3.SO3_wedge([0.1;-0.2;0.3]*2), [0.1;0.2;-0.3]*5;zeros(1,4)];
% W0 = [SE3_se3.SO3_wedge(rand(3,1)), rand(3,1);zeros(1,4)];
dt = 0.001;
num_steps = 1/dt*2;
% num_steps = 100;
J = diag([3,3,3,0.1,0.1,0.1]);
% J = diag([1,1,1,1,1,1]);
T = [0;0;0];
M = [0;0;0];
% T = rand(3,1)*4;
% M = rand(3,1)*5;



E1=E0;
W1=W0;
for ii = 1:num_steps



[Ed_1,Wd_1] = multirotor_dynamics(E1,W1,T,M,J);
[Ed_2,Wd_2] = multirotor_dynamics(E1+dt/2*Ed_1,W1+dt/2*Wd_1,T,M,J);
[Ed_3,Wd_3] = multirotor_dynamics(E1+dt/2*Ed_2,W1+dt/2*Wd_2,T,M,J);
[Ed_4,Wd_4] = multirotor_dynamics(E1+dt*Ed_3,W1+dt*Wd_3,T,M,J);

E1 = E1 + dt/6*(Ed_1 + 2*Ed_2 + 2*Ed_3 + Ed_4);
W1 = W1 + dt/6*(Wd_1+2*Wd_2 + 2*Wd_3 + Wd_4);

end

%% RK-MK

E2=E0;
W2=W0;
W2 = SE3_dse3.SE3_wedge(J*SE3_dse3.SE3_vee(W2));
Y = SE3_dse3(E2,W2);


for ii = 1:num_steps



y1 = dynamics2(Y,T,M,J,W0);
k1 = y1;

u2 = k1.*dt.*(0.5);
y2 = dynamics2(Y*u2.exp(),T,M,J,W0);
k2 = se3_dse3.dexp_inv(u2,y2);

u3 = k2.*dt.*(0.5);
y3 = dynamics2(Y*u3.exp(),T,M,J,W0);
k3 = se3_dse3.dexp_inv(u3,y3);

u4 = k3.*dt;
y4 = dynamics2(Y*u4.exp(),T,M,J,W0);
k4 = se3_dse3.dexp_inv(u4,y4);

v = (k1./6+k2./3+k3./3+k4./6).*dt;
% v = (y1./6+y2./3+y3./3+y4./6).*dt;
Y = Y*v.exp();

end

SE3_dse3.SE3_wedge(inv(J)*SE3_dse3.SE3_vee(Y.v));
TR = E0*expm(W0*2);

%% True values for constant translational acceleration