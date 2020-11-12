H = rand(3,3);
x = rand(2,1);
v = rand(2,1);


dt = 1e-5;
x1 = x;
x1(1) = x1(1) + dt;
x2 = x;
x2(2) = x2(2) + dt;

(constructG(H,x1)*v-constructG(H,x)*v)/dt
(constructG(H,x2)*v-constructG(H,x)*v)/dt

constructDG(H,x,v)

% H = sym('H',[3,3],'real');
% x = sym('x',[2,1],'real');
% v = sym('v',[2,1],'real');
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
