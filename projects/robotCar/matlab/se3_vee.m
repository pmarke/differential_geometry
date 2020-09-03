function V = se3_vee(M)
% Converts an element of se(3) to R^6
V=zeros(6,1);
w = M(1:3,1:3);
V(1:3) = so3_vee(w);
V(4:6) = M(1:3,4);

end

