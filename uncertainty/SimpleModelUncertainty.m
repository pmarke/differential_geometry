Q = 2;

x0=rand(4,1);

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];

dt = 0.1;
t1 = 0.5;
t2 = 2;
t3 = 4;
n = 1;
px1 = [];
xs1 = [];
px2 = [];
xs2 = [];
px3 = [];
xs3 = [];
for ii = 1:20000
  
x = 2*sqrt(Q*t1)*randn(1,1);
xs1 = [xs1, x];
p = (2*pi)^(-n/2)*det(Q*t1)^(-1/2)*exp( -1/2*x*inv(Q*t1)*x);
px1 = [px1,p ]  ;
    
end


for ii = 1:20000
  
x = 2*sqrt(Q*t1)*randn(1,1);
xs2 = [xs2, x];
p = (2*pi)^(-n/2)*det(Q*t2)^(-1/2)*exp( -1/2*x*inv(Q*t2)*x);
px2 = [px2,p ]  ;
    
end

for ii = 1:20000
  
x = 2*sqrt(Q*t1)*randn(1,1);
xs3 = [xs3, x];
p = (2*pi)^(-n/2)*det(Q*t3)^(-1/2)*exp( -1/2*x*inv(Q*t3)*x);
px3 = [px3,p ]  ;
    
end

figure(1), clf;
plot(xs1,px1,'*')
hold on
plot(xs2,px2,'*')
plot(xs3,px3,'*')
ylabel('p(x)')
legend('t=0.5','t=2','t=4')
xlabel('x')



%%
Q = 0.1;

xs = [];
for ii = 1:50000
 x = sqrt(Q)*randn(1,1);
 xs = [xs, x];   
end

figure(1),clf;
histogram(xs,100,'Normalization','pdf');