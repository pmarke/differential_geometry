




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

dt = 0.000001;

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

hh = [eye(3), hr1, hr2, hr3, hrho1, hrho2, hrho3, hw1, hw2, hw3]

% d1 = 0.1;
% d2 = 0.2;
% d3 = 0.3;
% 
h = H(R,w,p)
h - hh
% w = w*0
% 
% O = [H(R,d1*w,d1*p); H(R,d2*w,d2*p); H(R,d3*w,d3*p)]

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

% th^2*eye(3) +ssm(w)^2
% w*w'
% pV(w)*ssm(ssm(p)*w)
% myV(w+[0;0;1]*dt)
% myV(w)*myV(pV(w)*[0;0;1]*dt)
% myV(w) + pV(w)*dt



%%
syms a b c x y z phi real

w = [x;y;z];
p = [a;b;c];
th = sqrt(w'*w);


T = ssm(w)^3*p;

dT = [diff(T,x),diff(T,y),diff(T,z)];
dT = simplify(dT);
simplify(subs(dT,th^2,phi^2))

test = simplify(-2*ssm(w)^2*ssm(p)-ssm(ssm(w)^2*p)+ssm(w)*ssm(p)*ssm(w))

% ssm(ssm(w)^2*p)
% ssm(p)*ssm(w)^2
% ssm(w)*ssm(p)*ssm(w)

ssm(w)*ssm(p)*w;
ssm(p)*ssm(w)*w-ssm(w)*ssm(p)*w;
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

function V = pV(u)

th = norm(u);
U = ssm(u);

if abs(th) < 0.000001
    V = eye(3);
else
    V = eye(3)/2 + (sin(th)/th^3 - cos(th)/th^2)*U + ( (1-cos(th))/th^4 - sin(th)/th^3 + 1/(2*th^2))*U^2;
end

end

function h = H(R,w,p)

h = [eye(3) R*ssm(-myV(w)*p) R*myV(w) R*-pV(w)*ssm(-p)];

end