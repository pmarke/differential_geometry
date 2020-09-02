function V = so3_vee(M)

% Converts and element of so(3) to R^3

V = zeros(3,1);
V(1) = M(3,2);
V(2) = M(1,3);
V(3) = M(2,1);

end

