




rng('default');

p = rand(3,1)*0;
R = expm(ssm(rand(3,1)));
v = rand(3,1);
w = rand(3,1);
a = rand(3,1)*0;
g = [0;0;-9.81]*0;
h = 1;

Pf = p + h*v + h^2/2*g + h^2*R*A(h*w)*a
vf = v + h*g + h^2*R*A(h*w)*a
Rf = R*expm(ssm(w)*h);


g0 = [R p; 0 0 0 1];
u = R'*v;
int = 10000;
dt = h/int;


U = zeros(4,4);
U(1:3,1:3) = ssm(w);
U(1:3,4) = R'*u;
t = 0;
for ii = 1:int
    
    
    g0 = g0*myexp( dt*U);
    Rt = g0(1:3,1:3);
    du =  a;
    u = u + dt*du;
    U(1:3,4) = u;
    t = t+dt;

end


g0

W = U;
W(1:3,4) = rand(3,1);


%%
w = rand(3,1);
p = rand(3,1);
W = ssm(w);
P = ssm(p);
U = [W P; zeros(3,3) W]
u = [p;w];

norm(u)
norm(U)

%%
th = rand(1)-0.5
TH = [0 -th; th 0];

x = [0 -1; 1 0];
I = eye(2);

W = ( cos(th)-1 )/th*x + sin(th)/th*I;
W2 = (cos(norm(th)) -1)/norm(th)^2*TH + sin(norm(th))/norm(th)*I;
W-W2

%% Measurement Jacobian Test
rng('default');

dt = 0.0001;

R = expm(ssm(rand(3,1)));
w = rand(3,1);
t = rand(3,1);
p = rand(3,1);
delta = [1;0;0];

expm(ssm(delta*dt));

hr1 = (R*expm(ssm(delta*dt))*myV(w)*p - R*myV(w)*p)/dt;
delta = [0;1;0];
hr2 = (R*expm(ssm(delta*dt))*myV(w)*p - R*myV(w)*p)/dt;
delta = [0;0;1];
hr3 = (R*expm(ssm(delta*dt))*myV(w)*p - R*myV(w)*p)/dt;
hrho1 =  (R*myV(w)*(p+[1;0;0]*dt) - R*myV(w)*p)/dt;
hrho2 =  (R*myV(w)*(p+[0;1;0]*dt) - R*myV(w)*p)/dt;
hrho3 = (R*myV(w)*(p+[0;0;1]*dt) - R*myV(w)*p)/dt;
hw1 = (R*myV(w+[1;0;0]*dt)*(p) - R*myV(w)*p)/dt;
hw2 = (R*myV(w+[0;1;0]*dt)*(p) - R*myV(w)*p)/dt;
hw3 = (R*myV(w+[0;0;1]*dt)*(p) - R*myV(w)*p)/dt;

hh = [eye(3), hr1, hr2, hr3, hrho1, hrho2, hrho3, hw1, hw2, hw3];

d1 = 0.1;
d2 = 0.2;
d3 = 0.3;
d4 = 0.4;
% 
h = H(R,w,p);
h - hh;
% w = w*0
% 
O = [H(R,d1*w,d1*p); H(R,d2*w,d2*p); H(R,d3*w,d3*p);H(R,d4*w,d4*p)]
rank(O([1:7,10:12],[1:12]))
% expm(ssm(-t))*myV(w)
% myV(w)*expm(ssm(-t))

%%
dt = 0.00001;
y1 = (myV(w+[1;0;0]*dt)*p-myV(w)*p)/dt

th = norm(w);
W = ssm(w);
y2 = (1-cos(th))/th^2*ssm(-p) + (th-sin(th))/th^3*( -2*W*ssm(p)+ssm(p)*W) + (sin(th)*th - 2*(1-cos(th)))/th^4*ssm(-p)*w*w' + (-2*th + 3*sin(th)-cos(th)*th)/th^5*W*ssm(-p)*w*w';
% 
y2

y3 = (1-cos(th))/th^3*ssm(p) + (w'*p)*(4*cos(th) - 4 + th*sin(th))/(2*th^4)*ssm(w) + (th-sin(th))/th^3*(ssm(w)*ssm(p) + ssm(p)*ssm(w)) + (w'*p)*(3*sin(th)-th*cos(th)-2*th)/th^5*ssm(w)^2

% th^2*eye(3) +ssm(w)^2
% w*w'
% pV(w)*ssm(ssm(p)*w)
% myV(w+[0;0;1]*dt)
% myV(w)*myV(pV(w)*[0;0;1]*dt)
% myV(w) + pV(w)*dt

%%
% rng('default');
w = (rand(3,1)-0.5);
th = norm(w);
V = myV(w)
v = V(:,1)-[1;0;0];
t = v'*v
th;
t/th;

(1-cos(th))^2/th^4*(th^2-w(1)^2)+(th-sin(th))^2/th^6*(th^4-w(1)^2*th^2)
m = v(1);
-m*(2-2*cos(th)-sin(th)*th)/(th^2-sin(th)*th)-m
%%
syms a b c x y z phi real

w = [x;y;z];
p = [a;b;c];
th = sqrt(w'*w);

ssm(w)^3

% T = ssm(w)^3*p;
% 
% dT = [diff(T,x),diff(T,y),diff(T,z)];
% dT = simplify(dT);
% simplify(subs(dT,th^2,phi^2))
% 
% test = simplify(-2*ssm(w)^2*ssm(p)-ssm(ssm(w)^2*p)+ssm(w)*ssm(p)*ssm(w))

% ssm(ssm(w)^2*p)
% ssm(p)*ssm(w)^2
% ssm(w)*ssm(p)*ssm(w)

% ssm(w)*ssm(p)*w;
% ssm(p)*ssm(w)*w-ssm(w)*ssm(p)*w;
% T = w*w'*ssm(p)*ssm(w);
% simplify(T)

% ssm(ssm(-w)*p)-ssm(w)*ssm(p)
% 
% th = sqrt(w'*w);
% % 
% t1 = (1-cos(th))/th^2;
% pretty(simplify(diff(t1,x)))
% 
% t2 = (th - sin(th))/th^3;
% pretty(simplify(diff(t2,x)))
% 
% T = (1-cos(th))/th^2*ssm(w)*p + (th - sin(th))/th^3*ssm(w)^2*p+eye(3)*p;
% subs(simplify(diff(T,x)),th,phi)

%% Seeding test

% rng('default')
t = (rand(3,1)-0.5)*10;
phi = rand(1);
th = (rand(1)-0.5)*2;
psi = (rand(1)-0.5)*2;
R = rot(phi,th,psi);
w = (rand(3,1)-0.5)*2;
% w(3)=0;
% w(1)=0;
p = [(rand(1,1))*15;0;0];
% p = (rand(3,1)-0.5)*30;

g = [R t; zeros(1,3) 1];
u = [ssm(w) p; zeros(1,4)];

dt = 0.1;

g1 = g*expm(u*dt);
g2 = g*expm(u*2*dt);

y0 = t;
y1 = g1(1:3,4);
y2 = g2(1:3,4);


td0 = (y1-y0)/dt;
td1 = (y2-y1)/dt;

psi0 = atan2(td0(2),abs(td0(1)));
th0 = atan(-td0(3)*sin(psi0)/td0(2));
phi0 = 0;
R0 = rot(phi0,th0,psi0);
Rx = R0(:,1);
px = Rx'*td0/(Rx'*Rx)

% R*p
% p(1)
% 
% th
% th0
% psi
% psi0


tmp = R0'*td1/px/dt;
ye = -tmp(3);
ze = tmp(2);




% ye-w(2)
% ze - w(3)

w0 = [0;ye;ze]
% fun = @(xx) myfunc( R0'*td1/px ,dt,xx);
% fun = @(xx) myfunc( expm(dt*ssm(w)) ,dt,xx);
% options = optimoptions('fsolve','algorithm','levenberg-marquardt');
% ww = fsolve(fun,w0,options);
% ww
w

% tmp = R0'*td1/px;
% ze =sqrt( -2*(tmp(1)-1+W0(1,3)^2/2*dt^2)/dt^2);
% abs(ze)-abs(w(3))
% xe = 2*(tmp(2)-ze*dt)/W0(1,3)/dt^2
% xe = 2*(tmp(3)+ye)/ze/dt

% C = [1, -ye*3/(dt*ze), 6/dt^3/ze*(tmp(2)-ze*dt + (ze*ye^2+ze^3)*dt^3/6)];
% roots(C)

ge0 = [R0 y0; zeros(1,3) 1];
ue0 = [ssm(w0) [px;0;0]; zeros(1,4)];

ge0*expm(ue0*3)
g*expm(u*3)

function X = myexp(u)

w = u(1:3,1:3);
v = u(1:3,4);

th = norm(w);
R = eye(3) + sin(th)*w + (1-sin(th))*w^2;

if th < 0.000001
    B = eye(3);
else
    B = eye(3) + (1-cos(th))/th^2*w + (1-sin(th)/th)*w^2/th^2;
end

X = [R B*v; zeros(1,3) 1];

end


function J = Jl_inv(u)

th = norm(u);
U = wedge(u);

if th < 0.0000001
    J = eye(3)
else

J = eye(3) - 0.5*U + (2-th*cot(th/2))/(2*th^2)*U^2;
end
end



function u = A(v)

th = norm(v);
a = v/th;

x = (1-cos(th))/th^2*eye(3);
y = (th-sin(th))/th^2*ssm(a);
z = (th^2+2*cos(th)-2)/th^2*a*a';
u = x+y+z;

end


function W = ssm(w)

x = w(1);
y = w(2);
z = w(3);

W = [0 -z y;...
     z 0 -x;...
     -y x 0];


end


function V = myV(u)

th = norm(u);
U = ssm(u);

if abs(th) < 0.000001
    V = eye(3);
else
   V =  eye(3) + (1-cos(th))/th^2*U + (th-sin(th))/th^3*U^2;
end

end

function V = pV(w,p)

th = norm(w);
W = ssm(w);

if abs(th) < 0.000001
    V = eye(3)*ssm(p);
else
    V = (1-cos(th))/th^2*ssm(-p) + (th - sin(th))/th^3*(-2*W*ssm(p) + ssm(p)*W) + (sin(th)*th +2*(cos(th)-1))/th^4*W*p*w' + (3*sin(th)-cos(th)*th -2*th)/th^5*W^2*p*w';
end

end

function h = H(R,w,p)

h = [eye(3) R*ssm(-myV(w)*p) R*myV(w) R*pV(w,p)];

end

function R = rot(phi,th,psi)


R = [cos(th)*cos(psi)   sin(phi)*sin(th)*cos(psi)-cos(phi)*sin(psi)   cos(phi)*sin(th)*cos(psi)+sin(phi)*sin(psi);...
     cos(th)*sin(psi)   sin(phi)*sin(th)*sin(psi)+cos(phi)*cos(psi)   cos(phi)*sin(th)*sin(psi)-sin(phi)*cos(psi);...
     -sin(th)             sin(phi)*cos(th)                             cos(phi)*cos(th)];

end