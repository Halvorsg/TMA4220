syms('E','v');

C_inv = [1/E,-v/E,0;-v/E,1/E,0;0,0,(1+v)/E];
Ctemp = inv(C_inv);
C = Ctemp/(E/(1-v^2));
C = simplify(C);