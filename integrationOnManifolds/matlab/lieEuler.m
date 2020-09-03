
% Generate a random vector field on S2
v = rand(3,1);
V = so3_ssm(v);

fig = figure(1);
hold on
clf;

% Generate evenly spaced out points on the unit sphere in 3-d
theta=linspace(0,2*pi,20);
phi=linspace(0,pi,20);
[theta,phi]=meshgrid(theta,phi);
rho=1;
x=rho*sin(phi).*cos(theta);
y=rho*sin(phi).*sin(theta);
z=rho*cos(phi);
mesh(x,y,z);

scale = 0.2;


% Draw vector field all over
XX =x(:);
YY = y(:);
ZZ = z(:);

for ii=1:length(ZZ)
    p0 = [XX(ii),YY(ii),ZZ(ii)]';
    p1 = scale*V*p0+p0;
    mArrow3(p0,p1,'tipWidth',0.02,'stemWidth',0.01);
end

%% Draw flow from a single point
p0 = [1,0,0]';
pts1 = [];
for jj = 0:1:10

    pts1 = [pts1,expm(jj*V)*p0];
    
end
hold on;
plot3(pts1(1,:),pts1(2,:),pts1(3,:),'b');

p0 = [1,0,0]';
pts2 = [];
for jj = 0:0.01:10

    pts2 = [pts2,expm(jj*V)*p0];
    
end
hold on;
plot3(pts2(1,:),pts2(2,:),pts2(3,:),'r');

% mArrow3([0 0 0],[1 1 1],'color','red','stemWidth',0.02,'facealpha',0.5)