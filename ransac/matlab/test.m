%% Testing my model likelihood

rng('default')
PD = 0.9;
PG = 0.95;

lT = 0.9;
lF = 1.6;

lo = 0;

count = 0;

for kk = 1:100
ll = 0;
for jj = 1:10

    if rand(1) > PD*PG
        mT = 0;
    else
        mT = 1;
    end
    
%     mT = poissrnd(lT);
    mF = poissrnd(lF);
    
m = mT+mF;

%     s = 0;
%     for ii = 0:m
% 
%         s = s+ poisson(lT,ii)*poisson(lF, m-ii);
% 
%     end
%     s;
%     % ln = log(s/poisson(lF,m));
if m >0
    ln = log(PD*PG/lF*m+1-PD*PG);
else
    ln = log(1-PD*PG);
end

    ll = ln+ll;

end

ll;
p = 1 - 1/(1+exp(ll));
% p = 1 - 1/(1+10^ll);

if p < 0.5
    count = count + 1;
    
end

end

count


%%

% n = 4;
% PG = 0.8;
% g = chi2inv(PG,n);
% S = eye(3);
% c = pi^(3/2)/gamma(n/2+1);
% V = c*g^(n/2)*det(S)

%%
x = rand(1)*5;
h = rand(1)*2;
func((2*x+2*h))
funcT(2*x,h)

%%

syms x y z th real
assumeAlso(x^2+y^2+z^2 == th);

a = (1-cos(th))/th^2;
b = (th - sin(th))/th^3;

t =  [1;0;0]+a*[0 -z y]' + b*[-z^2-y^2 x*y x*z]';
simplify(t'*t)


function y = func(x)

y = 2*x^2 + x  +1;

end

function y = funcT(x,h)

y = func(x) + 2*(4*x +1)*h + 4*2*h^2;

end







function p = poisson(lambda, k)

p = exp(-lambda)*lambda^k/factorial(k);

end









