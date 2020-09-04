


syms w1 w2 v1 v2 u1 u2 real

W1 = [0 -w1;...
      w1 0];
  
W2 = [0 -w2;...
      w2 0];
  
  
r1 = [v1;v2];
r2 = [u1;u2];

X1 = [W1 r1; zeros(1,3)];
X2 = [W2 r2; zeros(1,3)];

X1*X2 - X2*X1

XX = [W1 [v2;-v1]; zeros(1,3)];
XX*[r2;w2]









%% Derivative test

rng('default')






function J = Jl_inv(v)

p = v(1:2);
th = v(3);

s = [0 -1; 1 0];

if abs(th) < 0.001
   J = eye(3); 
else
    
    W = (1-cos(th))/th*s + sin(th)/th*eye(2);
    D = (cos(th)-1)/th^2*s + (th-sin(th))/th^2*eye(2);
    
    W_i = inv(W);
    J = [W_i, -W_i*D*p;...
        zeros(1,2) 1];
    
end