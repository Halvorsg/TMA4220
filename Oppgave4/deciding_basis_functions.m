function phis = deciding_basis_functions()
syms('x','y','z')
A = [0,0,0,1;
     1,0,0,1;
     0,1,0,1;
     0,0,1,1];
 B = [x,y,z,1];
 vec = inv(A)*eye(4);
 phis = transpose(B*vec);
end