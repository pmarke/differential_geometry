syms w1 w2 w3 j11 j12 j13 j22 j23 j33 'real'
W = SE3_se3.SO3_wedge([w1;w2;w3])
J = [j11 j12 j13; j12 j22 j23; j13 j23 j33]
assumeAlso(J == J');
w = [w1;w2;w3];

J*w
