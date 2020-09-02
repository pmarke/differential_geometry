function V_wedge = se3_wedge(V)
% converts a vector V=(v;w) into the wedge form
% where v is translational and w is angular veloctiy.

V_wedge=[ssm(V(4:6)),V(1:3);zeros(1,4)];



end

