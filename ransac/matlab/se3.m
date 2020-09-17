




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