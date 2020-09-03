function w = ssm(v)

w1 = v(1);
w2 = v(2);
w3 = v(3);

w = [0 -w3 w2;...
    w3  0  -w1;...
    -w2 w1 0];

end