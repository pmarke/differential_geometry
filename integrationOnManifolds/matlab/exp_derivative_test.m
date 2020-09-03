v = rand(3,1); % translational velocity;
w = rand(3,1); % rotational velocity;

V = ssm(v);
W = ssm(w);

e = [v;w];
E = se3_wedge(e);


%% SO(3) tests

% for dexp positive and negative
sumn = 0;
sump = 0;
iter = 300;

for ii=0:iter
    
   tmp = W^ii/factorial(ii+1);
   sumn = sumn + (-1)^ii*tmp;
   sump = sump + tmp;
    
end

closedn = dexpwn(W);
closedp = dexpwp(W);
errorn = norm(sumn-closedn);
errorp = norm(sump-closedp);
max_error = 1e-9;
if (errorn > max_error)
    disp('so3 dexp n test failed')
end
if (errorp > max_error)
    disp('so3 dexp p test failed')
end

% for inverse dexp positive and negative 

if (norm(dexpwn_inv(W)*closedn-eye(3)) > max_error)
   disp('so3 inv dexp n test failed') 
end
if (norm(dexpwp_inv(W)*closedp-eye(3)) > max_error)
   disp('so3 inv dexp n test failed') 
end



%% SE3 tests

adx = se3_adx(e);

% for dexp positive and negative
sumn = 0;
sump = 0;

for ii=0:iter
    
   tmp = adx^ii/factorial(ii+1);
   sumn = sumn + (-1)^ii*tmp;
   sump = sump + tmp;
    
end

closedn = dexpEn(E);
closedp = dexpEp(E);
errorn = norm(sumn-closedn);
errorp = norm(sump-closedp);
max_error = 1e-9;
if (errorn > max_error)
    disp('se3 dexp n test failed')
    errorn
end
if (errorp > max_error)
    disp('se3 dexp p test failed')
    errorp
end

% for inverse dexp positive and negative 

if (norm(dexpEn_inv(E)*closedn-eye(6)) > max_error)
   disp('se3 inv dexp n test failed') 
   norm(dexpEn_inv(E)*closedn)
end
% if (norm(dexpwp_inv(W)*closedp-eye(3)) > max_error)
%    disp('so3 inv dexp n test failed') 
% end






%% Test


W1 = -ssm(rand(3,1))*2;
W2 = ssm(rand(3,1));
t = 0.1;
R1 = expm(W1);
R2 = expm(t*W2);
RT = R1*R2
% 
% expm(W1+dexpwn_inv(W1)*W2*t)

h = 1/10;
[q,r] = qr(magic(3));
Po = q;
Pv =Po;
Pe = Po;
error = [];

s = 0;
for ii=h:h:1

k1 = W1*Pv;
k2 = W1*(Pv+h/2*k1);
k3 = W1*(Pv+h/2*k2);
k4 = W1*(Pv+h*k3);
Pv = Pv + h/6*(k1+2*k2+2*k3+k4);





% O2 = F1/2;
% A2 = -ssm(h*W1*expm(O2)*Pe);
% F2 = dexpwp_inv(O2)*A2;
% 
% O3 = F2/2;
% A3 = -ssm(h*W1*expm(O3)*Pe);
% F3 = dexpwp_inv(O3)*A3;
% 
% O4 = F3;
% A4 = -ssm(h*W1*expm(O4)*Pe);
% F4 = dexpwp_inv(O4)*A4;

F1 = h*W1;

O2 = F1/2;
A2 = h*W1;
F2 = log_inv(O2,A2);

O3 = F2/2;
A3 = h*W1;
F3 = log_inv(O3,A3);

O4 = F3;
A4 = h*W1;
F4 = log_inv(O4,A4);

Ot = F1/6 + F2/3 + F3/3 + F4/6;
s = s+Ot;
Pe = expm(Ot)*Pe;

% k1 = h*W1*Pe;
% k2 = h*W1*expm(ssm(k1)/2)*Pe;
% k3 = h*W1*expm(ssm(k2)/2 - ssm(ssm(k1)*k2)/8)*Pe;
% k4 = h*W1*expm(ssm(k3))*Pe;
% 
% Pe = expm(1/6*(ssm(k1+2*k2+2*k3+k4-ssm(k1)*k4/2)))*Pe;
% tmp = expm(ii*W1)*Po;
% error = [error, tmp(1)-Pe(1)];
end
Pt = expm(W1)*Po
Pv
Pe
norm(Pt-Pe)

plot(error)

%% Another test

% W1 = ssm(rand(3,1))*10;
% W2 = ssm(rand(3,1));
t = 0.1;
R1 = expm(W1);
R2 = expm(t*W2);
RT = R1*R2
w1 = so3_vee(W1);
% 
% expm(W1+dexpwn_inv(W1)*W2*t)

h = 1/100;
[q,r] = qr(magic(3));
Po = expm(W2);
Pv =Po;
Pe = Po;
error = [];

Wo = logm(Pe);

s = 0;
for ii=h:h:1

% k1 = h*ssm(dexpwp_inv(Wo)*w1);
% k2 = h*ssm(dexpwp_inv(Wo+k1/2)*w1);
% k3 = h*ssm(dexpwp_inv(Wo+k2/2)*w1);
% k4 = h*ssm(dexpwp_inv(Wo+k3)*w1);

k1 = h*ssm(dexpwp_inv(Wo)*w1);
k2 = h*ssm(dexpwp_inv(logm(expm(k1/2)*expm(Wo)))*w1);
k3 = h*ssm(dexpwp_inv(logm(expm(k2/2)*expm(Wo)))*w1);
k4 = h*ssm(dexpwp_inv(logm(expm(k3)*expm(Wo)))*w1);



Ot = 1/6*(k1+2*k2+2*k3+k4);
s = s+Ot;
Wo = logm(expm(Ot)*expm(Wo));


end
Pt = expm(1*W1)*Po

Pe = expm(Wo)
norm(Pt-Pe)

plot(error)



%%
rng('default')
[q,r] = qr(magic(3));
yo = q;
yn = yo;
ye = yo;
yg = yo;
yf = yo;
yu = yo;
h = 1/5;

his = so3_vee(logm(q));
for ii = h:h:10
rng('default') 
k1 = vf(yn)*yn;
k2 = vf(yn+h/2*k1)*(yn+h/2*k1);
k3 = vf(yn+h/2*k2)*(yn+h/2*k2);
k4 = vf(yn+h*k3)*(yn+h*k3);
yn = yn + h/6*(k1+2*k2+2*k3+k4);  

    
rng('default')
u1 = h*vf(yu);

Y2 = u1/3;
W2 = vf(expm(Y2)*yu);
u2 = h*ssm(dexpwp_inv(Y2)*so3_vee(W2));
% u2 = h*W2;

Y3 = u2*2/3;
W3 = vf(expm(Y3)*yu);
u3 = h*ssm(dexpwp_inv(Y3)*so3_vee(W3));
% u3 = h*W3;

Y4 = u1/12 +u2/3 -u3/12;
W4 = vf(expm(Y4)*yu);
u4 = h*ssm(dexpwp_inv(Y4)*so3_vee(W4));
% u4 = h*W4;

Y5 = u1*25/48 -u2*55/24 + u3*35/48 + u4*15/8;
W5 = vf(expm(Y5)*yu);
u5 = h*ssm(dexpwp_inv(Y5)*so3_vee(W5));
% u5 = h*W5;

Y6 = u1*3/20 -u2*11/24 -u3/8 +u4/2 + u5/10;
W6 = vf(expm(Y6)*yu);
u6 = h*ssm(dexpwp_inv(Y6)*so3_vee(W6));
% u6 = h*W6;

Y7 = -u1*261/260 +u2*33/13 +u3*43/156 -u4*118/39 +u5*32/195 +u6*80/39;
W7 = vf(expm(Y7)*yu);
u7 = h*ssm(dexpwp_inv(Y7)*so3_vee(W7));
% u7 = h*W7;

YT = u1*13/200 +u3*11/40 + u4*11/40 +u5*4/25 + u6*4/25 +u7*13/200;
yu = expm(YT)*yu;

his = [his,so3_vee(logm(yu))];

rng('default')
F1 = h*vf(ye);  

O2 = F1/2;
A2 = vf(expm(O2)*ye);
F2 = h*ssm(dexpwp_inv(O2)*so3_vee(A2));
% F2 = h*log_inv(O2,A2);

O3 = F2/2;
A3 = vf(expm(O3)*ye);
F3 = h*ssm(dexpwp_inv(O3)*so3_vee(A3));% log_inv(O3,A3);
% F3 = h*log_inv(O3,A3);

O4 = F3;
A4 = vf(expm(O4)*ye);
F4 = h*ssm(dexpwp_inv(O4)*so3_vee(A4)); %loginv(O4,A4);
% F4 = h*log_inv(O4,A4);

Ot = F1/6 + F2/3 + F3/3 + F4/6;
ye = expm(Ot)*ye;

rng('default')
A1 = h*vf(yg);
B1 = A1;
A2 = h*vf(expm(B1/2)*yg);
B2 = A2-A1;
A3 = h*vf(expm(B1/2+B2/2-commute(B1,B2)/8)*yg);
B3 = A3-A2;
A4 = h*vf(expm(B1+B2+B3)*yg);
B4 = A4-2*A2+A1;
Og = B1+B2+B3/3 + B4/6 - commute(B1,B2)/6 - commute(B1,B4)/12;
yg = expm(Og)*yg;

rng('default')
t1 = h*vf(yf);
t2 = h*vf(expm(t1/2)*yf);
t3 = h*vf(expm(t2/2 - commute(t1,t2)/8)*yf);
t4 = h*vf(expm(t3)*yf);
Of = t1/6 + t2/3 + t3/3 + t4/6 - commute(t1,t4)/12;
yf = expm(Of)*yf;
  

end

yn;
ye;
yg;
yf;

en = norm(logm(TT'*yn))
eu = norm(logm(TT'*yu))
ee = norm(logm(TT'*ye))
ef = norm(logm(TT'*yf))
eg = norm(logm(TT'*yg))

figure(1), clf;
plot3(his(1,:),his(2,:),his(3,:));

% norm(BB-yf);

function a = commute(X,Y)
a = X*Y-Y*X;
end

function a = log_inv(A,C) 

tmp = (A*C-C*A);
a= C-tmp/2+(A*tmp-tmp*A)/12-ssm(A^4*so3_vee(C))/720;


end

function a = vf(y)
% a = (y-y')/2;
a = 3*sin((y-y'));
a = a*y;
a = (a-a');
% b = rand(3,3)*1;
% b = b-b';
% b = [0 1 -2; -1 0 1; 2 -1 0]*100;
% a = b;
end


function a = ath(th)

a = (cos(th)-1)/th^2;

end

function b = bth(th)

b = (th-sin(th))/th^3;

end

function c = cth(th)

c = -sin(th)/th^3 +2*( (1-cos(th))/th^4);

end

function d = dth(th)

d = -2/th^4 + 3/th^5*sin(th)-cos(th)/th^4;

end

function q = qn(th,w,v)

q = so3_vee(w)'*v*(cth(th)*w+dth(th)*w^2);

end

function q = qp(th,w,v)

q = so3_vee(w)'*v*(-cth(th)*w+dth(th)*w^2);

end

function d = dexpwn(w)

th = norm(w);

d = eye(3) + ath(th)*w+bth(th)*w^2;

end

function d = dexpwp(w)
th = norm(w);

d = eye(3) -ath(th)*w + bth(th)*w^2;

end

function d = dexpwn_inv(w)

th = norm(w);

d = eye(3) + w/2 - (th*cot(th/2)-2)/(2*th^2)*w^2;

end

function d = dexpwp_inv(w)
th = norm(w);

d = eye(3) - w/2 + w^2/12;

% d = eye(3) - w/2 - (th*cot(th/2)-2)/(2*th^2)*w^2;

end

function d = dexpEn(E)
w = E(1:3,1:3);
v = E(1:3,4);
th = norm(w);
dw = dexpwn(w);
q = so3_vee(w)'*v*(cth(th)*w+dth(th)*w^2);
b = ath(th)*ssm(v) + bth(th)*(w*ssm(v)+ssm(v)*w)+q;

d = [dw b; zeros(3,3) dw];

end

function d = dexpEp(E)
w = E(1:3,1:3);
v = E(1:3,4);
th = norm(w);
dw = dexpwp(w);
q = so3_vee(w)'*v*(-cth(th)*w+dth(th)*w^2);
b = -ath(th)*ssm(v) + bth(th)*(w*ssm(v)+ssm(v)*w)+q;

d = [dw b; zeros(3,3) dw];

end

function d = dexpEn_inv(E)
w = E(1:3,1:3);
v = E(1:3,4);
th = norm(w);
dw_inv = dexpwn_inv(w);
q = so3_vee(w)'*v*(cth(th)*w+dth(th)*w^2);
b = ath(th)*ssm(v) + bth(th)*(w*ssm(v)+ssm(v)*w)+q;

d = [dw_inv -dw_inv*b*dw_inv; zeros(3,3) dw_inv];

end

