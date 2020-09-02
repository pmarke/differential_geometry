function y = affineConnections(X,Y,gamma)
%UNTITLED14 Summary of this function goes here

syms x y z th phi Jw Jb mw mb l1 l2 mt 'real'

q = [x y th phi];

l = length(X);
y =sym(zeros(l,1));

for kk=1:l
    tmp = 0;
    for ii = 1:l
        tmp = tmp + diff(Y(kk),q(ii))*X(ii);
        for jj = 1:l
            tmp = tmp + gamma(kk,ii,jj)*X(ii)*Y(jj);
        end
    end 
    y(kk) = tmp;
end





end

