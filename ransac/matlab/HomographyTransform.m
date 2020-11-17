H = rand(3,3);
phi_h = rand(1,1);
t_h = rand(3,1);
n_h = [0;0;1];

R_h = [cos(phi_h) -sin(phi_h) 0; sin(phi_h) cos(phi_h) 0; 0 0 1];
H = R_h - t_h*n_h'

H1 = H(1:2,1:2);
h2 = H(1:2,3);
h3 = H(3,1:2);
h4 = H(3,3);
x = rand(2,1);
v = rand(2,1);
th = rand(1,1);
w = [0 -th; th 0];
o = [0 -1; 1 0];

B = constructB(H,x,v,w);
G = constructG(H,x);
wn = B*inv(G)

M = constructDG(H,x,v);
M*inv(G) + G*w*inv(G)


% t1 = M*v + G*w*v
% t2 = G*v

% 
% 
% 
% dt = 1e-5;
% x1 = x;
% x1(1) = x1(1) + dt;
% x2 = x;
% x2(2) = x2(2) + dt;
% G = constructG(H,x1);
% M = constructDG(H,x,v);
% % 
% k = [0 -1; 1 0];
% tmp = k*G*v;
% 
% th = tmp'*(M + G*w )*v/(tmp'*tmp)
% wn = [0 -th; th 0];
% 
% a = wn*G*v
% b = (M + G*w )*v

% (constructG(H,x1)*v-constructG(H,x)*v)/dt
% (constructG(H,x2)*v-constructG(H,x)*v)/dt

% constructDG(H,x,v)
% T = M + G*w
% wn = M*inv(G) + G*w*inv(G)
% wa = (wn-wn')/2;
% ws = (wn+wn')/2;

% H = sym('H',[3,3],'real');
% H1 = H(1:2,1:2);
% h2 = H(1:2,3);
% h3 = H(3,1:2);
% h4 = H(3,3);
% 
% x = sym('x',[2,1],'real');
% v = sym('v',[2,1],'real');
% syms th D 'real'
% w = [0 -th; th 0];
% D = (v'*G'*G*v)^1/2;
% D = norm(v'*G'*G*v);
% G = constructG(H,x);
% M = constructDG(H,x,v);
% e1 = [1;0];
% e2 = [0;1];
% 
% y = (M*v+G*w*v)/D - (G*(v*v')*G'*G*w*v)/D^3;
% z = (G*v)/D;
% tmp = e1'*( y);
% tmp2 = -e2'*z;
% tmp3 = simplify(tmp/tmp2)
% 
% tmp4 = simplify(-e2'*y/(-e1'*z))


% M_simp = constructDG_simp(H,x,v);

% v_new = constructG(H,x)*v;
% v_new = simplify(v_new)
% z = v_new(2)/v_new(1);
% z = simplify(z)
% 
% syms a b c 'real'
% assumeAlso(a == H(3,1:2)*x+H(3,3));
% 
% Blah = [diff(constructG(H,x)*v,x(1)), diff(constructG(H,x)*v,x(2))];
% subs(Blah,H(3,1:2)*x+H(3,3),a);
% Blah = simplify(Blah)
% 
% 
% 
% DG = constructDG(H,x,v);
% subs(DG,H(3,1:2)*x+H(3,3),a);
% DG = simplify(DG)
% 
% pretty(DG)


















function G = constructG(H,x)

H1 = H(1:2,1:2);
h2 = H(1:2,3);
h3 = H(3,1:2);
h4 = H(3,3);

tmp = h3*x+h4;
G = (tmp*H1 - (H1*x+h2)*h3)/tmp^2;

end

function DG = constructDG(H,x,v)
H1 = H(1:2,1:2);
h2 = H(1:2,3);
h3 = H(3,1:2);
h4 = H(3,3);
tmp = h3*x+h4;

DG = (2 * (H1*x + h2)*h3*v*h3 - (H1*(v*h3) + H1*(h3*v))*tmp)/tmp^3;

end

function DG = constructDG_simp(H,x,v)
H1 = H(1:2,1:2);
h2 = H(1:2,3);
h3 = H(3,1:2);
h4 = H(3,3);
tmp = h3*x+h4;
G = constructG(H,x);

DG = -G*(v*h3+eye(2)*(h3*v))/tmp;

end

function B = constructB(H,x,v,w) 
H1 = H(1:2,1:2);
h2 = H(1:2,3);
h3 = H(3,1:2);
h4 = H(3,3);
eta = h3*x+h4;


tmp1 = H1*w/eta
tmp2 = 2*(h3*v)*H1/eta^2
tmp3 = (H1*x+h2)*(h3*w)/eta^2
tmp4 = 2*(H1*x+h2)*(h3*v)*h3/eta^3

B = tmp1 - tmp2 - tmp3 + tmp4;


end

% function 
