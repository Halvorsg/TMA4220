function Ak = Stiffness_matrix()
syms('x_1','x_2','x_3','y_1','y_2','y_3','L1','L2','L3','x','y','v','real')

Lamb2Cart = [1,1,1;
            x_1,x_2,x_3;
            y_1,y_2,y_3];


        
        
Cart2Lamb = inv((Lamb2Cart));

K = [x_1-x_3,x_2-x_3;
    y_1-y_3,y_2-y_3];
A2 = det(K); A = det(K)/2;


compCart2Lamb = simplify(Cart2Lamb*A2);
disp(compCart2Lamb)

lambdas = compCart2Lamb*[1;x;y];
d  = [x,y];
D  = [d(1) , v;
      v    , d(2);
      d(2) , d(1)];
  
BasisFunctions = [eye(2)*lambdas(1) , eye(2)*lambdas(2) , eye(2)*lambdas(3)];
% B1 = [y_2-y_3 ,    0    , y_3-y_1 ,    0    , y_1-y_2 , 0;
%         0    , x_3-x_2 ,    0    , x_1-x_3 ,    0    , x_2-x_1;
%      x_3-x_2 , y_2-y_3 , x_1-x_3 , y_3-y_1 , x_2-x_1 , y_1-y_2];

B = diff(ones(3,6)*x,v);
for i = 1:3
    for j = 1:6
        for k = 1:2
            B(i,j) = B(i,j) + diff(BasisFunctions(k,j),D(i,k));
        end
    end
end
    
syms('E','v','x_','y_','real');
C_inv = [1 , -v , 0 ;
        -v , 1  , 0;
        0    ,     0,2*(1+v)]*1/E; %Be carefull with the last term here
Ctemp = inv(C_inv);
C = Ctemp/(E/(1-v^2));
C = simplify(C);
         
Element_Stiffness_Matrix = simplify(A*B'*C*B * E/(1-v^2)*1/(2*A)^2);

Ak = matlabFunction(Element_Stiffness_Matrix);


%% testing

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

         
         
         
         
         