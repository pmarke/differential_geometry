function Ad = Adm(X)

Ad = [X(1:2,1:2), -ssm(1)*X(1:2,3); 0 0 1];

end