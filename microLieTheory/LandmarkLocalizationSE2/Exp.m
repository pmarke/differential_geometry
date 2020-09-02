% This is the exponential map from the the cartesian vector 
% space to the SE(2) manifold.
% u: the input containing the angular velocity and liner velocity
% u = [vx,vy, w];

function M = Exp(u)

% Decompress the input u
vd = u(1:2);
wd = u(3);


% Construct the skew symmetric matrix with wd
ssm_wd = [0 -wd; wd 0];

% Construct the rotation matrix
R = [cos(wd) -sin(wd); sin(wd) cos(wd)];

% Construct the skew symmetric matrix with the identity element 1
ssm_i = [0 -1; 1 0];

% If norm(w) -delta_e > 0 for delta_e > 0, then 
% you can use the alternate method for the exponenetial map
% which is more computationally efficient. Else
% you will need to use the original method
if(norm(wd) > 0.00001)

    % Construct V
    V = sin(wd)*eye(2)/wd+ (1-cos(wd))/wd*ssm_i;

    % Construct the exponenetial map
    M = [R,  V*vd;...
         zeros(1,2)               1];
else
 
    m = [ssm_wd, vd; zeros(1,2) 0];
    M = expm(m);
end

end