function C = fun(x)

P1 = [1 0; 0 10e6];
P2 = [10e6 0; 0 1];

C = trace( inv(x*inv(P1) + (1-x)*inv(P2))     );




end