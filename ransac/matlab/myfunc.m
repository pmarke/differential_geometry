function e = myfunc(m,dt,v)
x = v(1);
y = v(2);
z = v(3);
th = norm(v*dt);

k = [1;0;0];
% n = [1;0;0] + [0;z;-y]*dt + [-z^2-y^2; x*y;x*z]*dt^2/2;

n = eye(3) + sin(th)/th*ssm(v*dt) + (1-cos(th))/th^2*ssm(v*dt)^2;


e = norm(m-n);


end


function W = ssm(w)

x = w(1);
y = w(2);
z = w(3);

W = [0 -z y;...
     z 0 -x;...
     -y x 0];


end