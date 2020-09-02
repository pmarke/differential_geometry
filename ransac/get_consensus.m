function   cons = get_consensus(x,measurements,t)

siz = size(measurements);

tau = 0.2;

cons = 0;

V = x(1:3);
P = x(4:6);
W = x(7:9);
z = [x(end);0;0];

RV = expm(ssm(V));

for jj = 1:length(t)
    
    yh = P -t(jj)*RV*Jr(W,-t(jj))*z;

    for ii = 1:siz(2)
        y = measurements(:,ii,length(t)-jj+1);
        e = (y-yh)'*(y-yh);
        if (e<tau)
            cons = cons+1;
        end
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