function M = so3_wedge(v)

% Takes a vector in R^3 can retruns the skew symmetric matrix of the vector

%Check if the input is a column vector 
if size(v,2)>1
v=v';
end

x = v(1);
y = v(2);
z = v(3);

M = [0 -z  y;...
     z  0 -x;...
    -y  x  0];

end

