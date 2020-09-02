




















o = observability(1,2,3)
rank(o)













function O = observability(alpha,w,x)

sa = sin(alpha);
ca = cos(alpha);

O = [ zeros(2,1)  eye(2)  zeros(2,1)  zeros(2,1);...
      -sa*x zeros(1,2) 0 ca;...
      ca*x zeros(1,2) 0 sa;...
      -ca*x*w zeros(1,2) -sa*x -sa*w;...
      sa*x*w zeros(1,2) -ca*x -ca*w;...
      sa*x*w^2 zeros(1,2) -2*ca*x*w -ca*w^2;...
      ca*x*w^2 zeros(1,2) w*sa*x*w sa*w^2];
    




end