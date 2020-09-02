syms A B 'real'
A = sym('A',[2 2],'real');
B = sym('B',[2 2],'real');

ath = [2 1 7]';
bth = [1 2 3]'*0.01;
aa = rand(3,1);
bb = rand(3,1);

% a = [ath;aa];
a = [aa;ath];
% b = [bth; bb];
b = [bb;bth];

A = [ssm(ath) aa; 0 0 0 1];
B = [ssm(bth) bb; 0 0 0 1];
A_ad = [ssm(ath) ssm(aa); zeros(3,3) ssm(ath)];
B_ad = [ssm(bth) ssm(bb); zeros(3,3) ssm(bth)];

A_ad*B_ad-B_ad*A_ad
t = A_ad*b;
T = [ssm(t(4:6)), ssm(t(1:3));zeros(3,3) ssm(t(4:6))]

%% SE2

t1 = [1 2 0.1]';
t2 = [2 -2 0.5]';

t1_ad = [0 -t1(3) t1(2);...
         t1(3) 0 -t1(1);...
         0 0 0];
t2_ad = [0 -t2(3) t2(2);...
         t2(3) 0 -t2(1);...
         0 0 0];
     
T1 = [0 -t1(3) t1(1);...
         t1(3) 0 t1(2);...
         0 0 0];
     
T2 = [0 -t2(3) t2(1);...
         t2(3) 0 t2(2);...
         0 0 0];
     
t1_ad*t2_ad - t2_ad*t1_ad
T1*T2-T2*T1
t1_ad*t2
         


% A = ssm(a);
% B = ssm(b);
% % 
% % X1 = my_exp(A)*my_exp(B);
% % X2 = my_exp(A+Jr_inv(A)*B);
% % X1 = simplify(X1)
% % X2 = simplify(X2)
% jr_i = JR(a);
% jl = JL(B);
% X1 = expm(A)*expm(B);
% X2 = expm(ssm(a+jr_i*b));
% X3 = expm(A+B);
% X4 = expm(inv(jl)*A+B);
% norm(X1 - X2)
% norm(X1-X3)
% norm(X1-X4)




function x = JR(A)

phi = norm(A);
a = A/phi;
x = phi/2*cot(phi/2)*eye(3)+(1-phi/2*cot(phi/2))*(a*a')+phi/2*ssm(a);

end

function x = JL(B)

x = 0;
k_f = 200;
for k=0:k_f
   x = x+1.0/factorial(k+1)*(B)^k
end

end

function x = Jr_inv(A)

x = 0;
k_f = 11;
B = [1,-1/2,1/6,0,-1/30,0,1/42,0,-1/30,0,5/66,0,-691/2730];

for k=0:k_f
    x = x +B(k+1)*(-A)^k/factorial(k)
end
end



function x = my_exp(A)

x = 0;
k_f = 11;

for k=0:k_f
    x = x +A^k/factorial(k);
end

end