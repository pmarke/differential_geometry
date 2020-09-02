



yo = 4;
yh = 1;
y1= yo+yh;

fy(yo);
t1 = fy(yo+yh);
t2 = fy(yo) + fy1(yo)*yh + fy2(yo)*yh^2/2 + fy3(yo)*yh^3/6 + fy4(yo)*yh^4/factorial(4) + fy5(yo)*yh^5/factorial(5) + fy6(yo)*yh^6/factorial(6);


xo = 1/2*yo^2;
xh = 1/2*(yo+yh)^2-xo;
x1 = xo+xh;
fx(xo)
fx(xo+xh);
fx(xo)+fx1(xo)*xh+fx2(xo)*xh^2/2 + fx3(xo)*xh^3/6;



fx(phi(yo))

phi(fx(yo))







function x = phi(y)

x = y^5/4;

end



function f = fx(x)

f = x^3/2;

end

function f = fx1(x)

f = x^2 + x;

end

function f = fx2(x)

f = 2*x + 1;

end

function f = fx3(x)
f = 2;
end



function f = fy(y)
f = y^6/24+y^4/8;
end

function f = fy1(y)
f = 1/4*y^5+1/2*y^3;
end

function f = fy2(y)
f = 5/4*y^4 + 3/2*y^2;
end

function f = fy3(y)
f = 5*y^3 + 3*y;
end

function f = fy4(y)
f = 15*y^2 + 3;
end

function f = fy5(y)
f = 30*y;
end

function f = fy6(y)
f = 30;
end