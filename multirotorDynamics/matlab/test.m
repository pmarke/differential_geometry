w = sym('w',[3,1],'real');
v = sym('v', [3,1],'real');

P = [ssm(w) ssm(v); zeros(3,3) ssm(w)];

I = rand(6,6)