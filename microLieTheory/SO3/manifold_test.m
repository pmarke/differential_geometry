w_x = 2*pi/20;
w_y = 2*pi/20;
w_z = 2*pi/20;
Ts = 0.01;

R = eye(3);
Rm = R;
Re = R;
for tt = Ts:Ts:2.5
ww = [0,0,w_z];  
wx = ssm(ww);
Re = Re + Re*wx*Ts;
Rm = Rm*expm(wx*Ts);
end

for tt = Ts:Ts:2.5
ww = [0,w_y,0];  
wx = ssm(ww);
Re = Re + Re*wx*Ts;
Rm = Rm*expm(wx*Ts);   
    
end

Rm

R_true = Rot_i_b(0,pi/4,pi/4)'

%%
w_x = 2*pi/20;
w_y = 2*pi/6;
w_z = 2*pi/5;
Ts = 0.01;
t_i = 1;
tf = t_i + Ts;

Td = derivativeTransfrom(w_x*t_i,w_y*t_i);
ww = inv(Td)*[w_x;w_y;w_z];

R_i = Rot_i_b(w_x*t_i, w_y*t_i,w_z*t_i)';
R_truth = Rot_i_b(w_x*tf, w_y*tf,w_z*tf)'
wx = ssm(ww);
% R_i = expm(wx*Ts);
Rd = R_i*wx;
R_k = R_i + Rd*Ts

Rm = R_i*expm(wx*Ts)





function [T]= derivativeTransfrom(phi,theta)

T = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);...
     0 cos(phi)            -sin(phi);...
     0 sin(phi)/cos(theta), cos(phi)/cos(theta)];

end