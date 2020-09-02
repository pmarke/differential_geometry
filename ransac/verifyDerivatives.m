
th = 0.1;

s1 = 0;

for ii = 0:100000
   s1 =  s1 + (2*ii+2)*(-1)^ii*th^(2*ii)/factorial(2*ii+3);
end

c1 = sin(th)/th^3 - cos(th)/th^2;

s2 = 0;

for ii = 0:100000
    s2 = s2 + (2*ii+3)*(-1)^ii*th^(2*ii)/factorial(2*ii+4);
end

c2 = (1-cos(th))/th^4 - sin(th)/th^3+1/(2*th^2);
%%
rng('default')
v = rand(3,3);
v = v-v';

th = norm(v);
w = v/th;

s3 = 0;

v = 0.5*v;
for ii = 0:1000
    s3 = s3+ v^ii*(ii+1)/factorial(ii+2);
    
end

c3 = eye(3)/2 + (sin(th)/th^3 - cos(th)/th^2)*v + ( (1-cos(th))/th^4 - sin(th)/th^3 + 1/(2*th^2))*v^2

v = 0.5*v;
c4 = eye(3)/2 + (sin(th)/th^3 - cos(th)/th^2)*v + ( (1-cos(th))/th^4 - sin(th)/th^3 + 1/(2*th^2))*v^2