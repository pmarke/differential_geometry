


PD = 0.9;
PG = 0.95;

lT = 0.9;
lF = 1.6;

lo = 0;

count = 0;

for kk = 1:100
ll = 0;
for jj = 1:100

    if rand(1) > PD*PG
        mT = 0;
    else
        mT = 1;
    end
    
%     mT = poissrnd(lT);
    mF = poissrnd(lF);
    
m = mT*0+mF;

    s = 0;
    for ii = 0:m

        s = s+ poisson(lT,ii)*poisson(lF, m-ii);

    end
    s;
    % ln = log(s/poisson(lF,m));
    ln = log(PD*PG/lF*m+1-PD*PG);

    ll = ln+ll;

end

ll;
p = 1 - 1/(1+exp(ll));

if p < 0.5
    count = count + 1;
    
end

end

count

function p = poisson(lambda, k)

p = exp(-lambda)*lambda^k/factorial(k);

end









