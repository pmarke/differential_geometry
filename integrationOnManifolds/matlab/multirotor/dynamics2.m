function y = dynamics2(Y,T,M,J,W0)

E = Y.g;
W = Y.v; % This is Jw

% a = W0;
a = SE3_dse3.SE3_wedge(J\SE3_dse3.SE3_vee(W));

R = E(1:3,1:3);
e3 = [0;0;1];
F = [T+9.81*R'*e3*0;M];

b = SE3_dse3.SE3_wedge(F);

y = se3_dse3(a,b);

end