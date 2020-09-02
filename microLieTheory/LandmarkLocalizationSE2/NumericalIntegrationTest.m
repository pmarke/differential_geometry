X = [eye(2) zeros(2,1); 0 0 1];   % Initial position
tau = [10 0 pi/2]';

% Manifold

dt = 0.1;
t = 0:dt:4;
fin = length(t);
X_history = zeros(3,3,length(t));
for ii=2:fin
    X = X*Exp(dt*tau);
    X_history(:,:,ii)=X;
end

P = reshape(X_history(1:2,3,:),2,[]);

% First order Euler Integration

x = zeros(3,1);
x_history = zeros(3,length(t));

for ii=2:fin
    theta = x(3);
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    v = R*tau(1:2);
    x_dot = [v;tau(3)];
    x = x+dt*x_dot;
    x_history(:,ii) = x;
end

% rk4
rk = zeros(3,1);
rk_history = zeros(3,length(t));

for ii=2:fin
    h = dt;
    k1 = h*f(rk,tau);
    k2 = h*f(rk+k1/2,tau);
    k3 = h*f(rk+k2/2,tau);
    k4 = h*f(rk+k3,tau);
    rk = rk+(k1+2*k2+2*k3+k4)/6;
    
    
    
    rk_history(:,ii) = rk;
end

%%

radius = tau(1)/tau(3);
xx = cos(0:pi/12:2*pi)*radius;
yy = sin(0:pi/12:2*pi)*radius+radius;

filepath = '/home/mark/projects/autonomousSystems/project/figures'

f1 = figure(1);
clf;
hold on
plot(xx,yy)
plot(P(1,:),P(2,:))
plot(x_history(1,:), x_history(2,:))
plot(rk_history(1,:), rk_history(2,:))
legend('True','Manifold','Euler 1st Order', 'RK4')
xlabel('x pos (m)')
ylabel('y pos (m)')
title('Numerical Integration Test')

set(f1,'Units','Inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f1,[filepath,'/NumericalIntegration'],'-dpdf','-r0')



function rk_dot = f(x,tau)
    theta = x(3);
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    v = R*tau(1:2);
    rk_dot = [v;tau(3)];
end
