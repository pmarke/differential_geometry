


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