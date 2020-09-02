function tau = Log(M)


% R = M(1:2,1:2); % extract the rotation matrix
% p = M(1:2,3);   % extract the position vector
% th = logm(R);   % Get Theta

tmp = logm(M);

tau = [tmp(1:2,3);tmp(2,1)];


end

