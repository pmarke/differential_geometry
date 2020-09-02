syms x y z dx dy dth dphi dv th phi Jw Jb mw mb l1 l2 mt bt br 'real'
assumeAlso(mt == mb+2*mw);

% Generalized velocity and generalized inertia tensor of the body frame
Vb = [0 0 th x y 0]';
Ib = diag([0,0,Jb,mb, mb, 0]);
ix = [1 0 0]';
iy = [0 1 0]';

% Kinetic Energy of the body frame
KEb = Vb'*Ib*Vb/2;

% Calculate metric tensor for the body frame
Gb = calculateMetricTensor(KEb);

% Generalized velocity and generalized inertia tensor of the left wheel
rhol = Vb(4:6) + l1*so3_wedge(Vb(1:3))*ix + l2*so3_wedge(Vb(1:3))*iy;
Vl = [0 0 phi rhol']';
Il = diag([0,0,Jw,mw, mw, 0]);

% Kinetic Energy of the body frame
KEl = Vl'*Il*Vl/2;

% Calculate metric tensor for the body frame
Gl = calculateMetricTensor(KEl);

% Generalized velocity and generalized inertia tensor of the right wheel
rhor = Vb(4:6) + l1*so3_wedge(Vb(1:3))*ix - l2*so3_wedge(Vb(1:3))*iy;
Vr = [0 0 phi rhor']';
Ir = diag([0,0,Jw,mw, mw, 0]);

% Kinetic Energy of the body frame
KEr = Vr'*Ir*Vr/2;

% Calculate metric tensor for the body frame
Gr = calculateMetricTensor(KEr);

% Calculate metric tensor for the entire system
G = Gb + Gl + Gr;
assumeAlso(mt == mb+2*mw);
syms alpha c1 'real'
assumeAlso(c1 == 2*mw*l1^2 + 2*mw*l2^2 + Jb);
G = simplify(G);

% gammaG = christoffelMetric(G);
gammaG = sym(zeros(4,4,4));

%% distribution D
syms c2 c3 'real'
q = [x y th phi];
X1 = [0;0;0;1];
X2 = [cos(th)*c2; sin(th)*c2; sin(phi); 0];
% X2 = [cos(th)*(cos(phi)*l1+l2*sin(phi)); sin(th)*(cos(phi)*l1+l2*sin(phi)); sin(phi); 0];

assumeAlso(c3 == simplify(innerG(X2,X2,G)));
%% Generalized christoffel symbols
gammaH = christoffelGeneralized([X1,X2],gammaG,G);
assumeAlso(c3 == simplify(innerG(X2,X2,G)));
gammaH = simplify(gammaH);
%% Forces and Torques
F = [cos(th); sin(th); 0;
    0];
Fd = [-bt*cos(th)*dv; -bt*sin(th)*dv; 0; -br*dphi];
T = [0;0;0;1];

% Compute vector forces
Gs = simplify(inv(G));
YF = Gs*F;
YFd = Gs*Fd;
YT = Gs*T;

Yf = innerG(YF,X1,G)/innerG(X1,X1,G) + innerG(YF,X2,G)/innerG(X2,X2,G);
Yf = simplify(Yf);
Yfd = innerG(YFd,X1,G)/innerG(X1,X1,G) + innerG(YFd,X2,G)/innerG(X2,X2,G);
Yfd = simplify(Yfd);
Yt = innerG(YT,X1,G)/innerG(X1,X1,G) + innerG(YT,X2,G)/innerG(X2,X2,G);
Yt = simplify(Yt);
