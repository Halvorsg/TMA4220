function [ output_args ] = Stiffness_matrix_3D( input_args )
%STIFFNESS_MATRIX_3D Summary of this function goes here
%   Detailed explanation goes here
%% Get B-matrix
syms('v','x','y','z','l_1','l_2','l_3','l_4','x_1','x_2','x_3','x_4','y_1','y_2','y_3','y_4','z_1','z_2','z_3','z_4','dx','dy','dz','real')

Lamb2Cart = [1  ,  1  ,  1  ,  1 ;
            x_1 , x_2 , x_3 , x_4;
            y_1 , y_2 , y_3 , y_4;
            z_1 , z_2 , z_3 , z_4];
Cart2Lamb = inv(Lamb2Cart);

K = [x_1 - x_4 , x_2 - x_4 , x_3 - x_4;
    y_1 - y_4 , y_2 - y_4 , y_3 - y_4;
    z_1 - z_4 , z_2 - z_4 , z_3 - z_4];
A6 = det(K);

compCart2Lamb = simplify(Cart2Lamb*A6);
disp(compCart2Lamb)

lambdas = compCart2Lamb*[1;x;y;z];
% for i = 1:4
%     for j = 1:3
%         der = simplify(diff(lambdas(i),xbar(j)));
%         fprintf('dLambda_%i/d%s = %s\n', i,char(xbar(j)),char(der))
%     end
% end

%e = D*u
d = [x,y,z];
D = [d(1)  ,  v    ,   v;
      v    , d(2)  ,   v;
      v    ,  v    ,  d(3);
      v    ,  d(3) ,  d(2);
      d(3) ,   v   ,  d(1);
      d(2) ,  d(1) ,   v];
%Using v instead of 0 such that we get 0 when we derivate by v;
BasisFunctions = [eye(3)*lambdas(1) , eye(3)*lambdas(2) , eye(3)*lambdas(3) , eye(3)*lambdas(4)];
%Making symbolic matrix of zeros
B = diff(ones(6,12)*x,v);
for i = 1:6
    for j = 1:12
        for k = 1:3
            B(i,j) = B(i,j) + diff(BasisFunctions(k,j),D(i,k));
        end
    end
end
%% Get C-matrix
syms('E','v','x_','y_','real');
C_inv = [1 , -v , -v , 0     , 0     , 0 ;
        -v , 1  , -v , 0     , 0     , 0 ;
        -v , -v ,  1 , 0     , 0     , 0 ;
        0  ,  0 ,  0 , 2+2*v , 0     , 0 ;
        0  ,  0 ,  0 , 0     , 2+2*v , 0 ;
        0  ,  0 ,  0 , 0     , 0     , 2+2*v]*1/E;
Ctemp = inv(C_inv);
C = Ctemp,%/(E/((1+v)*(1-2*v)));
C = simplify(C);
disp(C)

%% Element stiffness matrix

Element_Stiffness_Matrix = simplify((B'*C*B))*A6/6*(1/A6)^2;


end
