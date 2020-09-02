P = sym('P', [3,1],'real');
v = sym('v', [3,1], 'real');
w = sym('w', [3,1], 'real');
x= sym('x', [3,1], 'real');



zd = [Jr_inv(v)*w; expm(ssm(v))*x; zeros(6,1)];

F = myG(v,x)*Jr_inv(v)*w;
Fd = [diff(F,v(1)),diff(F,v(2)), diff(F,v(3)), diff(F,P(1)),diff(F,P(2)), diff(F,P(3)),diff(F,w(1)),diff(F,w(2)), diff(F,w(3)), diff(F,x(1)),diff(F,x(2)), diff(F,x(3))];
% Fd = simplify(Fd);

G = Fd*zd;
Gd = [diff(G,v(1)),diff(G,v(2)), diff(G,v(3)), diff(G,P(1)),diff(G,P(2)), diff(G,P(3)),diff(G,w(1)),diff(G,w(2)), diff(G,w(3)), diff(G,x(1)),diff(G,x(2)), diff(G,x(3))];


A = [P;...
    expm(ssm(v))*x;...
    F;...
    Fd*zd;...
    Gd*zd];

Ad = [diff(A,v(1)),diff(A,v(2)), diff(A,v(3)), diff(A,P(1)),diff(A,P(2)), diff(A,P(3)),diff(A,w(1)),diff(A,w(2)), diff(A,w(3)), diff(A,x(1)),diff(A,x(2)), diff(A,x(3))];



%%

p1 = 1;
p2 = 2;
p3 = 3;
v1 = 1;
v2 = 2;
v3 = 3;
w1 = 1;
w2 = 0;
w3 = 0;
x1 = 1;
x2 = 0;
x3 = 0;


y = double(subs(Ad))




rank(y(1:12,1:10))


%%
p1 = 1;
p2 = 2;
p3 = 3;
v1 = 0;
v2 = 0;
v3 = 0;
w1 = 0;
w2 = 1;
w3 = 1;
x1 = 2;
x2 = 0;
x3 = 0;

V = [v1;v2;v3];
w = [w1;w2;w3];
P = [p1;p2;p3];
x = [x1;x2;x3];
R = expm(ssm(V));
E = [R P; zeros(1,3) 1];
e = [ssm(w) x; zeros(1,4)];


t = 0:0.1:10;
a = zeros(3,length(t));
for ii = 1:length(t)
   a(:,ii) = t(ii)*R*Jr(ssm(w),t(ii))*x; 
end

Ta = P + a;

Ta(:,end)

figure(1),clf;
plot3(Ta(1,:),Ta(2,:),Ta(3,:))

axis('square')
xlabel('x');
ylabel('y');


%%
w = [0.1;0.2;0.3];
W = ssm(w);
t1 = 0.1;
t2 = 0.2;

% Jr(W,t2)-Jr(W,t1)
% Jr(W,t2-t1)

expm(t1*W)*Jr(W,t2-t1)
Jr(W,t2)

function G = myG(v,p)

g = [p(2)*v(2) + p(3)*v(3), p(2)*v(1) - 2*p(1)*v(2), p(3)*v(1) - 2*p(1)*v(3);...
     p(1)*v(2) - 2*p(2)*v(1),   p(1)*v(1) + p(3)*v(3), p(3)*v(2) - 2*p(2)*v(3);...
     p(1)*v(3) - 2*p(3)*v(1), p(2)*v(3) - 2*p(3)*v(2),   p(1)*v(1) + p(2)*v(2)];

th = sqrt(v'*v);
V = ssm(v);

G = sin(th)/th*ssm(-p) +(cos(th)/th^2 - sin(th)/th^3)*(V*p)*v' + (1-cos(th))/th^2*g + (sin(th)/th^3 - (1-cos(th))/th^4*2)*V^2*p*v';


end

function jr_inv = Jr_inv(v)

V = ssm(v);
th = norm(v);

jr_inv=th/2*cot(th/2)*eye(3) + (1-th/2*cot(th/2))*(v*v') + th/2*V;

end

function jr = Jr(w,t)
v = t*w;
th = norm(v);
jr = eye(3) + (1-cos(th))/th^2*v + (th - sin(th))/th^3*v^2;

end


function w = ssm(v)
x = v(1);
y = v(2);
z = v(3);

w = [0 -z y; z 0 -x; -y x 0];

end