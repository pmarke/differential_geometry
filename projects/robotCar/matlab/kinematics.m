syms dx dy dth dph th phi mb mw l1 l2 Jb Jw 'real'

pbd = [dx dy 0]';
prd = [-l1*sin(th)*dth + l2*cos(th)*dth;...
       + l1*cos(th)*dth + l2*sin(th)*dth;...
       0]+pbd;
   
pld = [dx-l1*sin(th)*dth - l2*cos(th)*dth;...
       dy+l1*cos(th)*dth - l2*sin(th)*dth;...
       0];
   
KE = (1/2)*mb*(pbd'*pbd) + (1/2)*mw*(prd'*prd) + (1/2)*mw*(pld'*pld);
KE = simplify(KE);

%% holonomic constraint
assumeAlso(-sin(th)*dx+cos(th)*dy==0);

c = -sin(th+phi)*(dx-l1*sin(th)*dth-l2*cos(th)*dth)+cos(th+phi)*(dy+l1*cos(th)*dth-l2*sin(th)*dth);
c = simplify(c)