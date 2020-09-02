function [Models,consensus] = ransac(measurements,dt,Q,R, iters,t,lambda)

index_time = 1:length(t);
index_meas = 1:lambda+3;
min_num = 3;

ti = flip(t);

A = [0,0,0,0,0,0,0,0,0,-1];


Models =zeros(10,lambda+3);
consensus = zeros(lambda+3,1);

% loop though each final measuremetn
for kk = lambda+1:lambda +3

best_cons = 0;
best_model = 0;
    
% try iter number of possibilites
for ii = 1:iters

% Select measurements
it = datasample(index_time(1:end-1), min_num);
it = [it,length(t)];
% im = randi(lambda+3,1,min_num);
im = randi(3,1,min_num)+lambda;
im = [im,kk];
% im = ones(1,min_num+1)*kk;

times = ti(it);
meas = zeros(3,min_num+1);
for jj = 1:min_num+1
    meas(:,jj) = measurements(:,im(jj),it(jj));
end
    

x0 = [randn(3,1)*2;meas(:,end);randn(3,1)*2;rand(1)];


fun = @(x)opt(x,times,meas,Q,R);
options = optimoptions('fmincon','MaxIterations',5000,'MaxFunctionEvaluations',6000);
[xf,fval] = fmincon(fun,x0,A,0,[],[],[],[],[],options);

cons = get_consensus(xf,measurements,t);

if (cons > best_cons)
   best_cons = cons;
   best_model = xf;
   fval
end

% if (fval < best_cons)
%    best_cons = fval;
%    best_model = xf;
% end

end

Models(:,kk) = best_model;
consensus(kk) = best_cons;

end
end