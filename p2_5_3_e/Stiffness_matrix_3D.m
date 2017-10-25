function [ output_args ] = Stiffness_matrix_3D( input_args )
%STIFFNESS_MATRIX_3D Summary of this function goes here
%   Detailed explanation goes here

syms('x','y','z','l_1','l_2','l_3','l_4','x_1','x_2','x_3','x_4','y_1','y_2','y_3','y_4','z_1','z_2','z_3','z_4','real')

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
xbar = [x;y;z];
for i = 1:4
    for j = 1:3
        der = simplify(diff(lambdas(i),xbar(j)));
        fprintf('dLambda_%i/d%s = %s\n', i,char(xbar(j)),char(der))
    end
end
        
        
    


end

