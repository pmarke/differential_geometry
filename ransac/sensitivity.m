




rng('default')

v = [1;2;-3];
v = ssm(v)
th = norm(v);
v = v/norm(v);

c = eye(3)/2 + (sin(th)/th^3 - cos(th)/th^2)*v + ( (1-cos(th))/th^4 - sin(th)/th^3 + 1/(2*th^2))*v^2

dt = 0.0000001;


p = rand(3,1);


d = (expm(v+ssm(dt*[1;0;0]))*p-expm(v)*p)/dt

c*ssm(p)
ssm(p)*c





x = 1;
y = 2;
z = 3;

th1 = sqrt(x^2 + y^2+z^2);
th2 = sqrt((x+dt)^2+y^2+z^2);
% 
% (sin(th2)/th2-sin(th1)/th1)/dt
% 
% cos(th1)/th1^2*x - sin(th1)/th1^3*x

((1-cos(th2))/(th2)^2 - (1-cos(th1))/th1^2)/dt

sin(th1)*x/th1^3 - (1-cos(th1))*2*x/th1^4


%% 
rng('default')
v = sym('v',[3,1],'real');
p = sym('p',[3,1],'real');
% th = sym('th',[1,1],'real');
% v = rand(3,1);
% p = rand(3,1);
th = sqrt(v'*v);
% th = sqrt(v'*v);
V = ssm(v);

G = [p(2)*v(2) + p(3)*v(3), p(2)*v(1) - 2*p(1)*v(2), p(3)*v(1) - 2*p(1)*v(3);...
     p(1)*v(2) - 2*p(2)*v(1),   p(1)*v(1) + p(3)*v(3), p(3)*v(2) - 2*p(2)*v(3);...
     p(1)*v(3) - 2*p(3)*v(1), p(2)*v(3) - 2*p(3)*v(2),   p(1)*v(1) + p(2)*v(2)];

d = sin(th)/th*ssm(-p) +(cos(th)/th^2 - sin(th)/th^3)*(V*p)*v' + (1-cos(th))/th^2*G + (sin(th)/th^3 - (1-cos(th))/th^4*2)*V^2*p*v'
% d = simplify(d);

% dt = 0.0000001;
% 
% 
% 
% 
% d2 = [(expm(V+ssm(dt*[1;0;0]))*p-expm(V)*p)/dt,(expm(V+ssm(dt*[0;1;0]))*p-expm(V)*p)/dt,(expm(V+ssm(dt*[0;0;1]))*p-expm(V)*p)/dt]

rank(sin(th)/th*ssm(-p) +(cos(th)/th^2 - sin(th)/th^3)*(V*p)*v' + (1-cos(th))/th^2*G + (sin(th)/th^3 - (1-cos(th))/th^4*2)*V^2*p*v')




g=th^2*v(1)*(-v(3)*p(2)+v(2)*p(3))-th^2*(v(2)*p(2)+v(3)*p(3)) +2*v(1)*( (-v(3)^2-v(2)^2)*p(1) + v(1)*v(2)*p(2) + v(1)*v(3)*p(3))

%%

rng('default')
v = sym('v',[3,1],'real');
p = sym('p',[1,1],'real');
p = [p;0;0];
th = sym('th',[1,1],'real');
% v = rand(3,1);
% p = [rand(1);0;0];
% th = sqrt(v'*v);
% th = sqrt(v'*v);
V = ssm(v);


G = [p(2)*v(2) + p(3)*v(3), p(2)*v(1) - 2*p(1)*v(2), p(3)*v(1) - 2*p(1)*v(3);...
     p(1)*v(2) - 2*p(2)*v(1),   p(1)*v(1) + p(3)*v(3), p(3)*v(2) - 2*p(2)*v(3);...
     p(1)*v(3) - 2*p(3)*v(1), p(2)*v(3) - 2*p(3)*v(2),   p(1)*v(1) + p(2)*v(2)];


d = sin(th)/th*ssm(-p) +(cos(th)/th^2 - sin(th)/th^3)*(V*p)*v' + (1-cos(th))/th^2*G + (sin(th)/th^3 - (1-cos(th))/th^4*2)*V^2*p*v';
d= simplify(d)






function w = ssm(v)
x = v(1);
y = v(2);
z = v(3);

w = [0 -z y; z 0 -x; -y x 0];

end