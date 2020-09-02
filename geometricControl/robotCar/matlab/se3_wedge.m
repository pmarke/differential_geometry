function M = se3_wedge(V)
% Maps vectors in R^6 to se3; (w,v)

%Check if the input is a column vector 
if size(V,2)>1
V=V';
end

w = so3_wedge(V(1:3));
M = [w, V(4:6); zeros(1,3) 1];
end

