syms a11 a12 a13 a21 a22 a23 a31 a32 a33 w1 w2 w3 real;



R = [a11 a12 a13; a21 a22 a23; a31 a32 a33];
W = [0 -w3 w2; w3 0 -w1; -w2 w1 0]; 
V = [w1;w2;w3];

% assumeAlso(det(R)==1);
assumeAlso(R*R' == eye(3))
assumeAlso(R'*R==eye(3));

T = simplify(R*W*R');




%%
Y = rand(3,3);
Y = Y-Y';
G = expm(Y);
h = rand(3,3);
h = h-h';
hh = [h(3,2);h(1,3);h(2,1)];
v = rand(3,1);
% G*h*G'
% 
% G*hh

G'*h*G
hh'*G