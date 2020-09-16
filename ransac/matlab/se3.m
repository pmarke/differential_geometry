




rng('default');

p = rand(3,1);
R = expm(ssm(rand(3,1)));
v = rand(3,1);
w = rand(3,1);
a = rand(3,1);
g = [0;0;-9.81];
h = 2;

Pf = p + h*v + h^2/2*g + h^2*R*A(h*w)*a
Rf = R*expm(ssm(w)*h);


g0 = [R p; 0 0 0 1];
u = [ssm(w) R'*v; zeros(1,4)];
int = 100;
dt = h/int;



for ii = 1:int
    

    g0 = g0*expm( dt*u);
    Rt = g0(1:3,1:3);
    du = [zeros(3,3) a+Rt'*g; zeros(1,4)];
    u = u + dt*du;
    
end


g0














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