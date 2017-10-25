function g = Stiffness_matrix_()
syms('x_2','x_1','x_3','y_1','y_2','y_3','L1','L2','L3','real')

Cart = [1,1,1;
    x_1,x_2,x_3;
    y_1,y_2,y_3];

lambdas = simplify(inv(Cart))*(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2);

K = [x_1-x_3,x_2-x_3;
    y_1-y_3,y_2-y_3];
A = det(K)/2;



B = [y_2-y_3 ,   0   , y_3-y_1 ,   0   , y_2-y_1 , 0;
               0   , x_3-x_2 ,   0   , x_1-x_3 ,   0   , x_2-x_1;
             x_3-x_2 , y_2-y_3 , x_1-x_3 , y_3-y_1 , x_2-x_1 , y_1-y_2];

B1 =1/(2*A)*[y_2-y_3 ,   0   , y_3-y_1 ,   0   , y_2-y_1 , 0;
               0   , x_3-x_2 ,   0   , x_1-x_3 ,   0   , x_2-x_1;
             1/2*(x_3-x_2) , 1/2*(y_2-y_3) , 1/2*(x_1-x_3) , 1/2*(y_3-y_1) , 1/2*(x_2-x_1) , 1/2*(y_1-y_2)];

    
syms('E','v','x_','y_','real');
C_inv = [1/E , -v/E , 0 ;
        -v/E , 1/E  , 0;
        0    ,     0,2*(1+v)/E];
Ctemp = inv(C_inv);
C = Ctemp/(E/(1-v^2));
C = simplify(C);
         
Element_Stiffness_Matrix = B'*C*B; % * E/(1-v^2)*1/(2*A)^2
g = matlabFunction(Element_Stiffness_Matrix);
         
         
         
         
         