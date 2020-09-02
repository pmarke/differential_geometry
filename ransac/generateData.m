%% Generate true targets

rng = ('default')

% Target 1 (Move in a circle)
V1 = [0;0;0]; % Attitude in the Lie algebra
W1 = [0;0.2;0.3]; % Angluar velocity in the body frame
P1 = [0;0;0]; % Position in the inertial frame
X1 = [1;0;0]; % Translational velocity in the body frame

% Target 2 (Move in a straight line)

V2 = [5;2;1]; % Attitude in the Lie algebra
W2 = [0.1;0;0]; % Angluar velocity in the body frame
P2 = [1;3;5]; % Position in the inertial frame
X2 = [2;0;0]; % Translational velocity in the body frame

% Target 3 (Move in a spiral)

V3 = [2;3;4]; % Attitude in the Lie algebra
W3 = [0.2;-0.2;0.1]; % Angluar velocity in the body frame
P3 = [2;2;2]; % Position in the inertial frame
X3 = [1.5;0;0]; % Translational velocity in the body frame

% Statistical Data
Q = 0.001*eye(3);    % Process noise covariance
q = chol(Q);
R = 0.001*eye(3);   % Measurement noise covariance
r = chol(R);
lambda = 10;        % Expected number of noise measurements per time step

% Other Parameters
dt = 0.1;           % Time step
t = 0:dt:10;        % Time array
region = 15;        % Survelliance region

% Generate Data
R1_array = zeros(3,3,length(t));
P1_array = zeros(3,length(t));
R2_array = R1_array;
P2_array = P1_array;
R3_array = R1_array;
P3_array = P1_array;

y1_array = P1_array;   % Measurements
y2_array = P2_array;
y3_array = P3_array;

noise_array = zeros(3,lambda,length(t));

for ii = 1:length(t)
    
    % generate noise
    w1 = q*randn(3,1);
    w2 = q*randn(3,1);
    w3 = q*randn(3,1);
    x1 = q*randn(3,1);
    x2 = q*randn(3,1);
    x3 = q*randn(3,1);
    
    v1 = r*randn(3,1);
    v2 = r*randn(3,1);
    v3 = r*randn(3,1);
%     w1 = q*randn(3,1)*0;
%     w2 = q*randn(3,1)*0;
%     w3 = q*randn(3,1)*0;
%     x1 = q*randn(3,1)*0;
%     x2 = q*randn(3,1)*0;
%     x3 = q*randn(3,1)*0;
%     
%     v1 = r*randn(3,1)*0;
%     v2 = r*randn(3,1)*0;
%     v3 = r*randn(3,1)*0;
    
    if (t(ii)==0)
        R1_array(:,:,ii) = expm(ssm(V1));
        P1_array(:,ii) = P1;
        y1_array(:,ii) = P1 + v1;
        
        R2_array(:,:,ii) = expm(ssm(V2));
        P2_array(:,ii) = P2;
        y2_array(:,ii) = P2 + v2;
        
        R3_array(:,:,ii) = expm(ssm(V3));
        P3_array(:,ii) = P3;
        y3_array(:,ii) = P3 + v3;
    else
        R1_array(:,:,ii) = R1_array(:,:,ii-1)*expm(dt*ssm(W1+w1));
        P1_array(:,ii) = P1_array(:,ii-1)+ dt*R1_array(:,:,ii-1)*Jr(W1+w1,dt)*(X1+x1);
        y1_array(:,ii) = P1_array(:,ii)+v1;
        
        R2_array(:,:,ii) = R2_array(:,:,ii-1)*expm(dt*ssm(W2+w2));
        P2_array(:,ii) = P2_array(:,ii-1)+ dt*R2_array(:,:,ii-1)*Jr(W2+w2,dt)*(X2+x2);
        y2_array(:,ii) = P2_array(:,ii)+v2;
        
        R3_array(:,:,ii) = R3_array(:,:,ii-1)*expm(dt*ssm(W3+w3));
        P3_array(:,ii) = P3_array(:,ii-1)+ dt*R3_array(:,:,ii-1)*Jr(W3+w3,dt)*(X3+x3);
        y3_array(:,ii) = P3_array(:,ii)+v3;
    end
    
    for jj = 1:lambda
       x = rand(1)*30 - 5;
       y = rand(1)*6 - 1;
       z = rand(1)*20 - 5;
       noise_array(:,jj,ii)=[x;y;z]; 
    end
    
    
end

measurements = zeros(3,lambda+3,length(t));
measurements(:,1:lambda,:) = noise_array;
measurements(:,lambda+1,:) = y1_array;
measurements(:,lambda+2,:) = y2_array;
measurements(:,lambda+3,:) = y3_array;

clusters = form_clusters(measurements,t,4,1);

%%
[Models,Cons] = ransac(measurements,dt,Q,R, 400,t,lambda)
% Rh = expm(ssm(M(1:3)));
% Ph = M(4:6);
% Wh = M(7:9);
% xh = M(end);
% 
% 
% Rf = R1_array(:,:,end);
% Pf = P1_array(:,end);
% 
% a = zeros(3,length(t));
% for ii = 1:length(t)
%    a(:,ii) = Ph - t(ii)*Rh*Jr(Wh,-t(ii))*[xh;0;0]; 
% end
% 
% % a = zeros(3,length(t));
% % for ii = 1:length(t)
% %    a(:,ii) = Pf - t(ii)*Rf*Jr(W1,-t(ii))*X1; 
% % end

%% Generate estimates
M1 = Models(:,lambda+1);
Rh1 = expm(ssm(M1(1:3)));
Ph1 = M1(4:6);
Wh1 = M1(7:9);
xh1 = M1(end);

a1 = zeros(3,length(t));
for ii = 1:length(t)
   a1(:,ii) = Ph1 - t(ii)*Rh1*Jr(Wh1,-t(ii))*[xh1;0;0]; 
end

M2 = Models(:,lambda+2);
Rh2 = expm(ssm(M2(1:3)));
Ph2 = M2(4:6);
Wh2 = M2(7:9);
xh2 = M2(end);

a2 = zeros(3,length(t));
for ii = 1:length(t)
   a2(:,ii) = Ph2 - t(ii)*Rh2*Jr(Wh2,-t(ii))*[xh2;0;0]; 
end

M3 = Models(:,lambda+3);
Rh3 = expm(ssm(M3(1:3)));
Ph3 = M3(4:6);
Wh3 = M3(7:9);
xh3 = M3(end);

a3 = zeros(3,length(t));
for ii = 1:length(t)
   a3(:,ii) = Ph3 - t(ii)*Rh3*Jr(Wh3,-t(ii))*[xh3;0;0]; 
end


%% Plot

figure(1),clf;

hold on;
plot3(y1_array(1,:),y1_array(2,:),y1_array(3,:),'g');
plot3(y2_array(1,:),y2_array(2,:),y2_array(3,:),'b');
plot3(y3_array(1,:),y3_array(2,:),y3_array(3,:),'k');

% plot3(a1(1,:),a1(2,:),a1(3,:),'g*')
% plot3(a2(1,:),a2(2,:),a2(3,:),'b*')
% plot3(a3(1,:),a3(2,:),a3(3,:),'k*')

for jj = 1:lambda
   plot3(squeeze(noise_array(1,jj,:)),  squeeze(noise_array(2,jj,:)),squeeze(noise_array(3,jj,:)), 'r*')
end

















function jr = Jr(w,t)
v = t*ssm(w);
th = norm(v);
jr = eye(3) + (1-cos(th))/th^2*v + (th - sin(th))/th^3*v^2;

end

function w = ssm(v)
x = v(1);
y = v(2);
z = v(3);

w = [0 -z y; z 0 -x; -y x 0];

end