% Create random elements of SE3
w1 = rand(3,3);
w1 = w1-w1';
R1 = expm(w1);

w2 = rand(3,3);
w2 = w2-w2';
R2 = expm(w2);

g1 = [R1, rand(3,1); zeros(1,3), 1];
g2 = [R2, rand(3,1); zeros(1,3), 1];
v1 = logm(g1);
v2 = logm(g2);

% Create random elements of se3
x1 = rand(3,3);
x1 = x1'-x1;
x2 = rand(3,3);
x2 = x2'-x2;
y1 = rand(3,1);
y2 = rand(3,1);

v1 = [x1, y1; zeros(1,4)];
v2 = [x2, y2; zeros(1,4)];

% Create SE3 XAd se3

H1 = SE3_se3(g1,v1);
H2 = SE3_se3(g2,v2);

% Test multiplication
H3 = H1*H2; 

% Test Inverse. There is machine precision error
H4 = H3*H3.inv();

% Test identity
H5 = H4.I();
H5.g;
H5.v;

% Test se3 matrix adjoint
a1 = v1*v2-v2*v1;
a2=SE3_se3.SE3_wedge(SE3_se3.se3_m_ad(v1)*SE3_se3.SE3_vee(v2));
norm(a1-a2);

%% Test the log function
[log_g,log_v] = H1.log()

B = [1, -0.5, 1/6, 0, -1/30, 0, 1/42, 7, -1/30, 0, 5/66, 0, -691/2730,0,7/6]; % Bernouli numbers
u = log_g;
w = H1.v;
s = 0;
k = length(B);

for ii = 0:k-1
  
  ad = w;
  for jj=1:ii
      ad = u*ad-ad*u;
  end
    
  s = s + (-1)^ii*B(ii+1)/factorial(ii)*ad;  
    
end
s

%% Test the exponential function

h1 = se3_se3(log_g,log_v);
Hk = h1.exp();
norm(inv(Hk.g)*H1.g)
norm(Hk.v-H1.v)

%% Test dexp_inv for se3_se3
[v,u] = H2.log();
h2 = se3_se3(v,u);

h3 = se3_se3.dexp_inv(h1,h2);

s_v = 0;
s_u = 0;
for ii = 0:k-1
  
  ad = h2;
  for jj=1:ii
      ad = se3_se3.lie_bracket(h1,ad);
  end
    
  s_v = s_v + (-1)^ii*B(ii+1)/factorial(ii)*ad.v;  
  s_u = s_u + (-1)^ii*B(ii+1)/factorial(ii)*ad.u;
    
end

norm(h3.v-s_v)
norm(h3.u-s_u)