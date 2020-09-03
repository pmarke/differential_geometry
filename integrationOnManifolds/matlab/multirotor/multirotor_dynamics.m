function [Ed,Wd] = multirotor_dynamics(E,W,T,M,J)


Ed = E*W;
R = E(1:3,1:3);
e3 = [0;0;1];
F = [T+9.81*R'*e3*0;M];


w_vee = SE3_dse3.SE3_vee(W);
w_ad = SE3_dse3.se3_m_ad(W);

Wd_vee = inv(J)*(F+ w_ad'*J*w_vee);
Wd = SE3_dse3.SE3_wedge(Wd_vee);


end

