function e = opt(x,t,meas,Q,R)

V = x(1:3);
P = x(4:6);
W = x(7:9);
z = [x(end);0;0];

RV = expm(ssm(V));
e = 0;

for ii = 1:length(t)
    y = meas(:,ii);
    yh = P -t(ii)*RV*Jr(W,-t(ii))*z;
%     Qb = Q_bar(W,t(ii),RV,Q)*0;

    if (t(ii) == 0)
         e = e + (y-yh)'*inv(R)*(y-yh)*100; 
    else
         e = e + (y-yh)'*inv(R)*(y-yh); 
    end
    
     
    
    
end






end

function jr = Jr(w,t)

v = t*ssm(w);
th = norm(v);
if (th < 0.001)
    jr = eye(3);
else
    jr = eye(3) + (1-cos(th))/th^2*v + (th - sin(th))/th^3*v^2;
end

end

function w = ssm(v)
x = v(1);
y = v(2);
z = v(3);

w = [0 -z y; z 0 -x; -y x 0];

end

function Qb = Q_bar(W,t,R,Q)

V = Jr(W,-t);
Qb = t^2*R*V*Q*V'*R';

end