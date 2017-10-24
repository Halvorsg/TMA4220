function Ak = Stiffness_matrix()
syms('x_1','x_2','x_3','y_1','y_2','y_3','L1','L2','L3','real')

% Cart = [1,1,1;
%     x_1,x_2,x_3;
%     y_1,y_2,y_3];

% lambdas = simplify(inv(Cart))*(x_1*y_2 - x_2*y_1 - x_1*y_3 + x_3*y_1 + x_2*y_3 - x_3*y_2);

K = [x_1-x_3,x_2-x_3;
    y_1-y_3,y_2-y_3];
A = det(K)/2;



B = [y_2-y_3 ,    0    , y_3-y_1 ,    0    , y_1-y_2 , 0;
        0    , x_3-x_2 ,    0    , x_1-x_3 ,    0    , x_2-x_1;
     x_3-x_2 , y_2-y_3 , x_1-x_3 , y_3-y_1 , x_2-x_1 , y_1-y_2]/(2*A);
 
% B1 =1/(2*A)*[y_2-y_3 ,   0   , y_3-y_1 ,   0   , y_2-y_1 , 0;
%                0   , x_3-x_2 ,   0   , x_1-x_3 ,   0   , x_2-x_1;
%              1/2*(x_3-x_2) , 1/2*(y_2-y_3) , 1/2*(x_1-x_3) , 1/2*(y_3-y_1) , 1/2*(x_2-x_1) , 1/2*(y_1-y_2)];

    
syms('E','v','x_','y_','real');
C_inv = [1/E , -v/E , 0 ;
        -v/E , 1/E  , 0;
        0    ,     0,(1+v)/E];
Ctemp = inv(C_inv);
C = Ctemp/(E/(1-v^2));
C = simplify(C);
         
Element_Stiffness_Matrix = simplify(A*B'*C*B * E/(1-v^2)*1/(2*A)^2);

Ak = matlabFunction(Element_Stiffness_Matrix);

testC =   8*[4,1,0;
           1,4,0;
           0,0,2];
testElement_Stiffness_Matrix = A*B'*testC*B;

Aktest = matlabFunction(testElement_Stiffness_Matrix);

Ke = [6,3,-4,-2,-2,-1;
      3,6,2,4,-5,-10;
      -4,2,24,-12,-20,10;
      -2,4,-12,24,14,-28;
      -2,-5,-20,14,22,-9;
      -1,-10,10,-28,-9,38;];
  disp(Aktest(0,3,2,0,1,2)-Ke);
% See https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch15.d/IFEM.Ch15.Slides.d/IFEM.Ch15.Slides.pdf
%  for derivation
end

         
         
         
         
         