
rng('default');

dt = 0.000001;

R = expm(ssm(rand(1)));
w = rand(1);
t = rand(2,1);
p = rand(2,1);
delta = 1;

expm(ssm(delta*dt))

% (R*expm(ssm(delta*dt))*myV(w)*p - R*myV(w)*p)/dt
% (R*myV(w)*(p+[1;0]*dt) - R*myV(w)*p)/dt
% (R*myV(w)*(p+[0;1]*dt) - R*myV(w)*p)/dt
% (R*myV(w+delta*dt)*(p) - R*myV(w)*p)/dt


d1 = 0.1;
d2 = 0.2;
d3 = 0.3;

h = H(R,w,p);
w = w*0

O = [H(R,d1*w,d1*p); H(R,d2*w,d2*p); H(R,d3*w,d3*p)]



rank(O(1:6,1:6))
rank(O([1:4,6],[1:4,6]))

%%
x = 0.1;
th = 0.2;

-x*cos(-th)
-x*sin(-th)



function g =  myExp(u)

p = u(1:2);
th = u(3);

R = [cos(th) -sin(th); sin(th) cos(th)];
V = myV(u);

g = [R V*p; zeros(1,3)];


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

h = [eye(2) R*ssm(1)*myV(w)*p R*myV(w) R*pV(w)*p];

end