













%% RKMK 4
rng('default')
h = 2;
w0 = rand(3,1);
g0 = eye(3);
dw = rand(3,1)/10;

u1 = zeros(3,1);
m1 = zeros(3,1);
k1 = w0;
t1 = dw;
k1t = k1;
t1t = t1;

u2 = h/2*k1t;
m2 = h/2*t1t;
[k2,t2] = vectorfield(g0,w0,u2,m2,dw);

[k2t,t2t] = JL_INV(u2,m2,k2,t2);

u3 = h/2*k2t;
m3 = h/2*t2t;
[k3,t3] = vectorfield(g0,w0,u3,m3,dw);
[k3t,t3t]=JL_INV(u3,m3,k3,t3);

u4 = h*k3t;
m4 = h*t3t;
[k4,t4] = vectorfield(g0,w0,u4,m4,dw);
[k4t,t4t]=JL_INV(u4,m4,k4,t4);

v1 = h*(k1t/6 + k2t/3 + k3t/3 + k4t/6);
v2 = h*(t1t/6 + t2t/3 + t3t/3 + t4t/6);

[g1,g2] = myexp(v1,v2);
[g1,g2] = mult(g0,w0,g1,g2);
gt = expm(wedge(w0*h) + wedge(dw)*h^2/2);
% gt = expm(wedge(w0*h) + Jl(w0*h)*dw);
% wt = w0 + h*dw;
% 
error = vee(logm(inv(g1)*gt))
% g2-wt;
%%
i = 500;
% h = 1;
% rng('default');
% g0 = eye(3)
% w0 = rand(3,1);
% dw = rand(3,1);
dt = h/i;
R = g0;
u = w0;
v = dw;
t = 0;
for ii = 1:i
%     u = u+dt*v;
%    Rd = R*wedge(u);
%    R = R+Rd*dt;
    R = expm(wedge(dt*u))*R;
   u = u+dt*v;
   t = t+dt;
end

vee(logm(inv(R)*g1))
vee(logm(inv(R)*gt))

% T = expm(h*wedge(u))*g0
% vee(logm(inv(R)*T))
% wt-u;
% t;
% Jl(w0*h)*dw
% R
% gt = expm(wedge(w0*h) + wedge(Jl(w0*h/2)*h/2*dw))

%%
% u
% Jl_inv(0.2*u)*u

%%
% inv(expm(wedge(u)))*Jr(u)
% Jr(u)*Jl_inv(u)*Jr(u)

function [k,t] = vectorfield(g0,w0,u,m,dw)

[g1,g2] = myexp(u,m);
[g1,g2] = mult(g0,w0,g1,g2);

k = g2;
t = dw;

end




function [g1,g2] = mult(g0,w0,n,m)

g1 = n*g0;
g2 = vee(n*wedge(w0)*inv(n)) +m;

% g1 = n*g0;
% g2 = w0 +m;

end


function [g1,g2] = myexp(u,v)

g1 = expm(wedge(u));
g2 = Jl(u)*v;

% g1 = expm(wedge(u));
% g2 = v;

end



function W = wedge(w)

x = w(1);
y = w(2);
z = w(3);

W = [0 -z y;...
     z 0 -x;...
     -y x 0];

end

function w = vee(W)

w = [W(3,2);W(1,3);W(2,1)];

end

function J = Jr(u)

U = wedge(u);
th = norm(u);

J = eye(3) + (cos(th)-1)/th^2*U + (th-sin(th))/th^3*U^2;

end

function J = Jl(u)

U = wedge(u);
th = norm(u);

J = eye(3) + (-cos(th)+1)/th^2*U + (th-sin(th))/th^3*U^2;

end

function J = Jr_inv(u)
U = wedge(u);
th = norm(u);

J = eye(3) + 0.5*U - (th*cot(th/2)-2)/(2*th^2)*U^2;

end


function J = JJl_inv(u,v)

th = norm(u);
U = wedge(u);
V = wedge(v);

a = -0.5*U;
b = (2-th*cot(th/2))/(2*th^2)*(U*V+V*U);
c = -(u'*v)*(8-th^2-2*th*cot(th/2)-th^2*cot(th/2)^2)/(4*th^4)*U^2;
J = a +b +c;

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

function [a,b] = JL_INV(u,v,x,y)

% a = Jl_inv(u)*x;
% b = JJl_inv(u,v)*x+Jl_inv(u)*y;

a = Jl_inv(u)*x;
b = y;

end
