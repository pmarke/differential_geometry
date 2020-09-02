rng('default')

H = [1/5,1/10,1/100,1/1000,1/10000];
det_his = zeros(1,length(H));
time_rk4_hist = zeros(1,length(H));
time_rkmk6_hist = zeros(1,length(H));
time_rkmk4_hist = zeros(1,length(H));
time_rkmk4opt_hist = zeros(1,length(H));

error_rk4_hist = zeros(1,length(H));
error_rkmk5_hist = zeros(1,length(H));
error_rkmk4_hist = zeros(1,length(H));
error_rkmk4opt = zeros(1,length(H));

for jj = 1:length(H)

[q,r] = qr(magic(3));
yo = q;
yn = yo;
ye = yo;
yf = yo;
yu = yo;
h = H(jj);

time_rk4 = 0;
time_rkmk6 = 0;
time_rkmk4 = 0;
time_rkmk4opt = 0;

his = so3_vee(logm(q));
for ii = h:h:10

tic;
    
% RK 4
k1 = vf(yn)*yn;
k2 = vf(yn+h/2*k1)*(yn+h/2*k1);
k3 = vf(yn+h/2*k2)*(yn+h/2*k2);
k4 = vf(yn+h*k3)*(yn+h*k3);
yn = yn + h/6*(k1+2*k2+2*k3+k4);  


time_rk4 = time_rk4 + toc;
  
tic;

% RK-MK 6
u1 = h*vf(yu);

Y2 = u1/3;
W2 = vf(my_exp(Y2)*yu);
u2 = h*ssm(dexpwp_inv(Y2)*so3_vee(W2));
% u2 = h*W2;

Y3 = u2*2/3;
W3 = vf(my_exp(Y3)*yu);
u3 = h*ssm(dexpwp_inv(Y3)*so3_vee(W3));
% u3 = h*W3;

Y4 = u1/12 +u2/3 -u3/12;
W4 = vf(my_exp(Y4)*yu);
u4 = h*ssm(dexpwp_inv(Y4)*so3_vee(W4));
% u4 = h*W4;

Y5 = u1*25/48 -u2*55/24 + u3*35/48 + u4*15/8;
W5 = vf(my_exp(Y5)*yu);
u5 = h*ssm(dexpwp_inv(Y5)*so3_vee(W5));
% u5 = h*W5;

Y6 = u1*3/20 -u2*11/24 -u3/8 +u4/2 + u5/10;
W6 = vf(my_exp(Y6)*yu);
u6 = h*ssm(dexpwp_inv(Y6)*so3_vee(W6));
% u6 = h*W6;

Y7 = -u1*261/260 +u2*33/13 +u3*43/156 -u4*118/39 +u5*32/195 +u6*80/39;
W7 = vf(my_exp(Y7)*yu);
u7 = h*ssm(dexpwp_inv(Y7)*so3_vee(W7));
% u7 = h*W7;

YT = u1*13/200 +u3*11/40 + u4*11/40 +u5*4/25 + u6*4/25 +u7*13/200;
yu = my_exp(YT)*yu;


time_rkmk6 = time_rkmk6 +toc;

% his = [his,so3_vee(logm(yu))];

tic;

% RK-MK 4
F1 = h*vf(ye);  

O2 = F1/2;
A2 = vf(my_exp(O2)*ye);
F2 = h*ssm(dexpwp_inv(O2)*so3_vee(A2));
% F2 = h*log_inv(O2,A2);

O3 = F2/2;
A3 = vf(my_exp(O3)*ye);
F3 = h*ssm(dexpwp_inv(O3)*so3_vee(A3));% log_inv(O3,A3);
% F3 = h*log_inv(O3,A3);

O4 = F3;
A4 = vf(my_exp(O4)*ye);
F4 = h*ssm(dexpwp_inv(O4)*so3_vee(A4)); %loginv(O4,A4);
% F4 = h*log_inv(O4,A4);

Ot = F1/6 + F2/3 + F3/3 + F4/6;
ye = my_exp(Ot)*ye;


time_rkmk4 = time_rkmk4 + toc;

tic;
% optimized
t1 = h*vf(yf);
t2 = h*vf(my_exp(t1/2)*yf);
t3 = h*vf(my_exp(t2/2 - commute(t1,t2)/8)*yf);
t4 = h*vf(my_exp(t3)*yf);
Of = t1/6 + t2/3 + t3/3 + t4/6 - commute(t1,t4)/12;
yf = my_exp(Of)*yf;


time_rkmk4opt = time_rkmk4opt + toc;
  

end

% time_rk4 
% time_rkmk6
% time_rkmk4 
% time_rkmk4opt 

det_his(jj) = det(yn);

time_rk4_hist(jj) = time_rk4;
time_rkmk6_hist(jj) = time_rkmk6;
time_rkmk4_hist(jj) = time_rkmk4;
time_rkmk4opt_hist(jj) = time_rkmk4opt;



% 
en = norm(logm(TT'*yn));
eu = norm(logm(TT'*yu));
ee = norm(logm(TT'*ye));
ef = norm(logm(TT'*yf));

error_rk4_hist(jj) = en;
error_rkmk5_hist(jj) = eu;
error_rkmk4_hist(jj) = ee;
error_rkmk4opt(jj) = ef;

% figure(1), clf;
% plot3(his(1,:),his(2,:),his(3,:));

jj
end
%%
figure(1),clf;
loglog(H,error_rkmk5_hist)
hold on
loglog(H,error_rk4_hist)
loglog(H,error_rkmk4_hist)
loglog(H,error_rkmk4opt)
legend('RK-MK 6','RK 4', 'RK-MK 4', 'RK-MK 4 Opt');

figure(2),clf;
loglog(H,time_rkmk6_hist)
hold on
loglog(H,time_rk4_hist)
loglog(H,time_rkmk4_hist)
loglog(H,time_rkmk4opt_hist)
legend('RK-MK 6','RK 4', 'RK-MK 4', 'RK-MK 4 Opt');

f3 = figure(3);
clf;
loglog(time_rkmk6_hist,error_rkmk5_hist,'-*')
hold on
loglog(time_rk4_hist,error_rk4_hist,'-*')
loglog(time_rkmk4_hist,error_rkmk4_hist,'-*')
loglog(time_rkmk4opt_hist,error_rkmk4opt,'-*')
legend('RK-MK 6','RK 4', 'RK-MK 4', 'RK-MK 4 Opt');
xlabel('Computation Time (s)')
ylabel('Error')

f4 = figure(4);
clf;
loglog(H,abs(det_his-1),'*-')
hold on
xlabel('Time Steps')
ylabel('Error')

% save("dataFast.mat","time_rk4_hist","time_rkmk6_hist","time_rkmk4_hist","time_rkmk4opt_hist","error_rk4_hist","error_rkmk5_hist","error_rkmk4_hist","error_rkmk4opt","det_his")

filepath = '/home/mark/thinkspace/researchnotes/geometry/integrationOnManifolds/images';

set(f3,'Units','Inches');
pos = get(f3,'Position');
set(f3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f3,[filepath,'/example2Error'],'-dpdf','-r0')

set(f4,'Units','Inches');
pos = get(f4,'Position');
set(f4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f4,[filepath,'/example2Det'],'-dpdf','-r0')

function a = commute(X,Y)
a = X*Y-Y*X;
end

function a = log_inv(A,C) 

tmp = (A*C-C*A);
a= C-tmp/2+(A*tmp-tmp*A)/12-ssm(A^4*so3_vee(C))/720;


end

function R = my_exp(W)

th = norm(W);
w = W/th;

R = eye(3) + sin(th)*w + (1-cos(th))*w^2;

end

function R = my_exp2(W)

th = norm(W);
w = W/th;

R = eye(3) + w + w^2/2 + w^3/6 + w^4/24;

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
