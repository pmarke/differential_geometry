
%% Measurement Jacobian Test
rng('default');

dt = 0.000001;

R = expm(ssm(rand(1)));
w = rand(1);
t = rand(2,1);
p = rand(2,1);
delta = 1;



u = [p;w];

getH(R,w)*getF(delta,u)


% (R*expm(ssm(delta*dt))*myV(w)*p - R*myV(w)*p)/dt
% (R*myV(w)*(p+[1;0]*dt) - R*myV(w)*p)/dt
% (R*myV(w)*(p+[0;1]*dt) - R*myV(w)*p)/dt
% (R*myV(w+delta*dt)*(p) - R*myV(w)*p)/dt


d1 = 0.1;
d2 = 0.2;
d3 = 0.3;

h = H(R,w,p)
% w = w*0
% 
% O = [H(R,d1*w,d1*p); H(R,d2*w,d2*p); H(R,d3*w,d3*p)]
% 
% 
% 
% rank(O(1:6,1:6))
% rank(O([1:4,6],[1:4,6]))



%%
% v = [1;2;0.2];
% jl = Jl(v);
% wl = Wl(v(3));
% dl = Dl(v(3));
% jr = Jr(v);
% 
% vee(wedge(v));
% 
% jl_inv = Jl_inv(v)
% Jr_inv(-v)
% jr_inv = Jr_inv(v);
% 
% dt = 1e-4;
% e1 = [1;0;0];
% e2 = [0;1;0];
% e3 = [0;0;1];
% tmp1 = inv(myExp(v))*expm(wedge(e1*dt));
% tmp2 = inv(myExp(v))*myExp(e2*dt);
% tmp3 = inv(myExp(v))*myExp(e3*dt);

% ( vee(logm(myExp(v)*myExp(e1*dt))) - v)/dt
% 
% ( vee(logm(myExp(v)*myExp(e2*dt))) - v)/dt
% 
% ( vee(logm(myExp(v)*myExp(e3*dt))) - v)/dt

function h = getH(R,w)
W = [0 -w; w 0];
    h = [R*expm(W) zeros(2,4)];
end

function F = getF(delta,  u) 
w = [0 -1; 1 0];

U = expm(wedge(-delta*u));
Ad = [U(1:2,1:2), w*U(1:2,3); zeros(1,2) 1];

F = [ Ad Jr(delta*u)*delta; zeros(3,3) eye(3)];

end

function W = wedge(u)

W = [ssm(u(3)), u(1:2); zeros(1,3)];

end

function w = vee(W)

w = [W(1:2,3);W(2,1)];

end

function J = Jl(u)
p = u(1:2);
th = u(3);

J = [Wl(th), Dl(th)*p; 0, 0,1];

end

function J = Jr(u)
p = u(1:2);
th = u(3);

J = [Wr(th), Dr(th)*p; 0, 0,1];

end

function J = Jr_inv(u)

p = u(1:2);
th = u(3);
W = inv(Wr(th));
J = [W, -W*Dr(th)*p; 0,0,1];

end

function J = Jl_inv(u)

p = u(1:2);
th = u(3);
W = inv(Wl(th));
J = [W, -W*Dl(th)*p; 0,0,1];

end


function g =  myExp(u)

p = u(1:2);
th = u(3);

R = [cos(th) -sin(th); sin(th) cos(th)];
V = Wl(th);

g = [R V*p; zeros(1,2),1];


end

function V = myV(u)

th = u;

if abs(th) < 0.000001
    V = eye(2);
else
   V =  sin(th)*eye(2)/th + (1-cos(th))/th*ssm(1);
end

end

function U = ssm(u)

U = [0 -u; u 0];

end

function V = pV(u)

th = u;

if abs(th) < 0.000001
    V = eye(2);
else
    V = ( cos(th)*th - sin(th) )/th^2*eye(2) + ( sin(th)*th - (1-cos(th)))/th^2*ssm(1);
end

end

function h = H(R,w,p)

h = [R R*ssm(1)*myV(w)*p R*myV(w) R*pV(w)*p];

end

function W = Wr(th)

if (abs(th) > 1e-7 )

a = (cos(th) - 1)/th;
b = sin(th)/th;

W = a*ssm(1) + b*eye(2);

else
W = eye(2);
end

end


function W = Wl(th)

if (abs(th) > 1e-7 )

a = (1-cos(th))/th;
b = sin(th)/th;

W = a*ssm(1) + b*eye(2);

else
W = eye(2);
end

end

function D = Dr(th)

if(abs(th) > 1e-7)
    
a = (1-cos(th))/th^2;
b = (th-sin(th))/th^2;

D = a*ssm(1) + b*eye(2);
    
else
   D = eye(2); 
end

end

function D = Dl(th)

if(abs(th) > 1e-7)
    
a = (cos(th)-1)/th^2;
b = (th-sin(th))/th^2;

D = a*ssm(1) + b*eye(2);
    
else
   D = eye(2); 
end

end