
rng('default')

Q = 0.1*eye(3);
R = 0.001*eye(3);
P = eye(6);

u = rand(3,1);
g = eye(3);
dt = 0.1;

v = rand(3,1)*4;
h = expm(wedge(rand(3,1)));



for ii = 1:100
   
    
    [P,g,u] = propagate(P,g,u,Q,dt);
    h = h*expm(wedge(dt*v));
    
    
    [P,g,u] = update(P,g,u,h,R);
    
end







function [P,g,u] = propagate(P,g,u,Q,dt)

Adm = Ad(expm(wedge(-u*dt)));
J = Jr(u*dt);

F = [Adm,J;zeros(3,3),eye(3,3)];

G = [J; zeros(3,3)];

P = F*P*F' +G*Q*G';

g = g*expm(wedge(u*dt));

end


function [P,g,u] = update(P,g,u,z,R)

e = vee(logm(inv(g)*z));

H = [eye(3), zeros(3,3)];

Z = H*P*H' + R;
K = P*H'*inv(Z);

delta = K*e;

g = g*expm(wedge(delta(1:3)));
u = u + 10*delta(4:6);
P
P = P-K*Z*K'

end



function w = ssm(x)

w = [0 -x; x 0];


end


function u = wedge(v)

    p = v(1:2);
    th = v(3);
    
    u = [ ssm(th) p; zeros(1,3)];

end

function v = vee(u)

p = u(1:2,3);
th = u(2,1);

v = [p;th];

end


function Adm = Ad(g)

R = g(1:2,1:2);
t = g(1:2,3);

s = [0 1; -1 0];

Adm = [R s*t; zeros(1,2), 1];

end





function J = Jl_inv(v)

p = v(1:2);
th = v(3);

if abs(th) < 0.001
   J = eye(3); 
else
    
    Q = ssm
    
end

end

function J = Jr(v)

v1 = v(1);
v2 = v(2);
th = v(3);

if abs(th) < 0.001
    J = eye(3);
else

J = [sin(th)/th (1-cos(th))/th (th*v1-v2+v2*cos(th)-v1*sin(th))/th^2;...
      (cos(th)-1)/th sin(th)/th (v1 +th*v2 -v1*cos(th) -v2*sin(th))/th^2;...
      0 0 1];
end

end