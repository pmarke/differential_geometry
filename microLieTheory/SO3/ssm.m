function M = ssm(V)
%SSM: converts a vector into a skew symmetric matrix
M=[0     -V(3)    V(2) ;...
    V(3)   0     -V(1) ;...
   -V(2)  V(1)    0 ];
end

