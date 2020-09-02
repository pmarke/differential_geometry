Q = diag([0.1,0.1,0.1]);    % Process noise covariance
R = diag([0.1,0.1]);         % Measurment noise covariance
P = eye(3);                  % Error covariance
P_history = [sqrt([P(1,1);P(2,2);P(3,3)])];

v = [0.002; 0; 0.01];          % Velocity

X = [eye(2) zeros(2,1); 0 0 1]; % Initial state


for ii = 1:70
    
    X = X*Exp(v);
    
    F = Adm(Exp(-v));
    G = Jr(v);
    
    P = F*P*F' + G*Q*G'*0;
    
    P_history = [P_history,sqrt([P(1,1);P(2,2);P(3,3)])];
    
end

v = -v;

for ii = 1:70
    
    X = X*Exp(v);
    
    F = Adm(Exp(-v));
    G = Jr(v);
    
    P = F*P*F' + G*Q*G'*0;
    
    P_history = [P_history,sqrt([P(1,1);P(2,2);P(3,3)])];
    
end

C = chol(P,'lower');
C = C(:,1:2);
CC = [C,-C,C(:,1)];
p = [];
for jj = 1:5
   
    delta = CC(:,jj);
    Dx = Exp(delta);
    t = X*Dx;
    p = [p,t(1:2,3)];
    
end


figure(1),clf
plot(p(1,:),p(2,:),'k')
hold on
plot(X(1,3),X(2,3),'b*')
xlim([-10,10])
ylim([-10,10])

figure(2), clf;

subplot(3,1,1);
plot(P_history(1,:))

subplot(3,1,2)
plot(P_history(2,:))
subplot(3,1,3)
plot(P_history(3,:))
