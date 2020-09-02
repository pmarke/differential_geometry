 syms a1 a2 p11 p12 p21 p22 'real'
 
 B1 = [0 -a1 p11; a1 0 p12; 0 0 0];
 B2 = [0 -a2 p21; a2 0 p22; 0 0 0];
 
 1/2*(B1*B2-B2*B1)
 