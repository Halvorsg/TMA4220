function TestCase()
syms('v','x','y','u','v','real');
%% 2D Case
% The function u
U = [(x^2-1)*(y^2-1);
    (x^2-1)*(y^2-1)];

Uxx = diff(diff(U(1),x),x);
Uxy = diff(diff(U(1),x),y);
Uyy = diff(diff(U(1),y),y);
Vxx = diff(diff(U(2),x),x);
Vxy = diff(diff(U(2),x),y);
Vyy = diff(diff(U(2),y),y);

gradS = [Uxx + (v+1)/2*Vxy + (1-v)/2*Uyy;
        Vyy + (v+1)/2*Uxy + (1-v)/2*Vxx];
gradS = simplify(gradS);
%% 2D matrix form
syms('v','real')
d = [x,y]; %Being used indirectly as the gradient
%D is the matrix of derivatives s.t. Du = e
D = [d(1)  ,  v    
      v    , d(2) 
      d(1) ,  d(2)];
% C is s.t. sigma = C*e
C =     [ 1, v,         0;
         v, 1,         0;
         0, 0, 1/2 - v/2];
% Pre allocating an empty syms matrix
B = diff(ones(3,1)*x,v);
% Finding the e's
for i = 1:3
    for j = 1
        for k = 1:2
            B(i,j) = B(i,j) + diff(U(k,j),D(i,k));
        end
    end
end
% Finding sigma on vector form
S = C*B;
% sigma on matrix form
SM = [S(1),S(3);S(3),S(2)];

%The gradient of sigma
dS = [diff(SM(1,1),x) + diff(SM(2,1),y);
      diff(SM(1,2),x) + diff(SM(2,2),y)];
% The gradient of sigma found in a for-loop
ddS = diff(dS,u);
for i = 1:2
    ddS(i) = diff(SM(1,i),x) + diff(SM(2,i),y);
end
disp(expand(dS))

%% 3D case
syms('z','real')
% The function u
U = [(x^2-1)*(y^2-1)*(z^2-1);
    (x^2-1)*(y^2-1)*(z^2-1);
    (x^2-1)*(y^2-1)*(z^2-1)];

% Creating the C matrix s.t. sigma = C*e
syms('E','v','x_','y_','real');
C_inv = [1 , -v , -v , 0     , 0     , 0 ;
        -v , 1  , -v , 0     , 0     , 0 ;
        -v , -v ,  1 , 0     , 0     , 0 ;
        0  ,  0 ,  0 , 2+2*v , 0     , 0 ;
        0  ,  0 ,  0 , 0     , 2+2*v , 0 ;
        0  ,  0 ,  0 , 0     , 0     , 2+2*v]*1/E;
Ctemp = inv(C_inv);
C = Ctemp/(E/((1+v)*(1-2*v)));
C = simplify(C);

% Creating indirect gradient and matrix of derivative operators
syms('v','real')
d = [x,y,z];
D = [d(1)  ,  v    ,   v;
      v    , d(2)  ,   v;
      v    ,  v    ,  d(3);
      v    ,  d(3) ,  d(2);
      d(3) ,   v   ,  d(1);
      d(2) ,  d(1) ,   v];
%Using v instead of 0 such that we get 0 when we derivate by v;
%Making symbolic matrix of zeros
B = diff(ones(6,1)*x,v);
% Making B, a vector containing e's
for i = 1:6
    for k = 1:3
        B(i,1) = B(i,1) + diff(U(k,1),D(i,k));
    end
end
% Finding sigma as a vector
S = C*B;

%Sigma in matrix form
SM = [S(1),S(6),S(5);
      S(6),S(2),S(4);
      S(5),S(4),S(3) ];
%Gradient of sigma
dS = [diff(SM(1,1),x) + diff(SM(2,1),y) + diff(SM(3,1),z);
      diff(SM(1,2),x) + diff(SM(2,2),y) + diff(SM(3,2),z);
      diff(SM(1,3),x) + diff(SM(2,3),y) + diff(SM(3,3),z)];

% Finding the gradient of sigma by a for-loop.
ddS = diff(dS,E);
for i = 1:3
    ddS(i) = diff(SM(i,1),x) + diff(SM(i,2),y) + diff(SM(i,3),z);
end
  
disp(simplify(ddS))

