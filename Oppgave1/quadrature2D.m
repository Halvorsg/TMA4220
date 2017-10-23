function I = quadrature2D(p1,p2,p3,Nq,g)
%Pre-conditions
    % p1 , p1 , p3 - Vectors defining the vertices of the triangle
    % Nq - Number of points evaluated
    % g(x,y) - function evaluated 
% Post-conditions
    % I - a vector of 3 values. One for each basisfunction g(x,y) is
    % muliplied with
    
p = @(lambda) lambda(1)*p1 + lambda(2)*p2 + lambda(3)*p3;
Area = 1/2*abs(det([p1-p3,p2-p3]));
I = 0;
[lambdas, weights] = GaussWeights2D(Nq);
for i = 1:Nq
    point = p(lambdas(i,:));
    I = I + weights(i)*g(point(1),point(2))*lambdas(i,:)';%lambdas(i,:) here defines the basis functions
end
I = Area*I;
end

function [lambdas, weights] = GaussWeights2D(Nq)
    switch Nq
        case 1
            lambdas = [1/3,1/3,1/3];
            weights = 1;
        case 3
            lambdas = [1/2 , 1/2 ,   0;
                       1/2 ,   0 , 1/2;
                       0   , 1/2 , 1/2];
                   
            weights = [1/3 , 1/3 , 1/3];
        case 4 
            lambdas = [1/3 , 1/3 , 1/3;
                       3/5 , 1/5 , 1/5;
                       1/5 , 3/5 , 1/5;
                       1/5 , 1/5 , 3/5];
            weights = [-9/16 , 25/48 , 25/48 , 25/48];
        case 6
            x1 = 0.81684757298; x2 = 0.09157621351; x3 = 0.09157621351;
            x4 = 0.10810301817; x5 = 0.44594849092; x6 = 0.44594849092;
            w1 = 0.10995174366; w2 = 0.22338158968;
            lambdas = [x1 , x2 , x3;
                       x3 , x1 , x2;
                       x2 , x3 , x1;
                       x4 , x5 , x6;
                       x6 , x4 , x5;
                       x5 , x6 , x4];
           weights = [w1 , w1 , w1 , w2 , w2 , w2];
        otherwise
            error('No such Nq at this moment\n choose between 1,3,4,6')            
    end
end
            
    