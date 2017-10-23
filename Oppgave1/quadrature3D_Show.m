function I = quadrature3D_Show(p1,p2,p3,p4,Nq,g)
%Pre-conditions
    % p1 , p1 , p3 , p4 - Vectors defining the vertices of the tetrahedra
    % Nq - Number of points evaluated
    % g(x,y,z) - function evaluated 
% Post-conditions
    % I - a vector of 4 values. One for each basisfunction g(x,y,z) is
    % muliplied with

p = @(lambda) lambda(1)*p1 + lambda(2)*p2 + lambda(3)*p3 + lambda(4)*p4;
Area = 1/6*abs(det([p1-p4,p2-p4,p3-p4]));
I = 0;
[lambdas, weights] = GaussWeights3D(Nq);
for i = 1:Nq
    point = p(lambdas(i,:));
    I = I + weights(i)*g(point(1),point(2),point(3));
end
I = Area*I;
end


function [lambdas, weights] = GaussWeights3D(Nq)
    switch Nq
        case 1
            lambdas = [0.25 , 0.25 , 0.25 , 0.25];
            weights = 1;
        case 4
            lambdas = [0.5854102, 0.1381966, 0.1381966, 0.1381966;
                       0.1381966, 0.5854102, 0.1381966, 0.1381966;
                       0.1381966, 0.1381966, 0.5854102, 0.1381966;
                       0.1381966, 0.1381966, 0.1381966, 0.5854102];
            weights = [0.25 , 0.25 , 0.25 , 0.25];
        case 5
            lambdas = [0.25 , 0.25 , 0.25 , 0.25;
                       0.5  , 1/6  , 1/6  , 1/6;
                       1/6  , 0.5  , 1/6  , 1/6;
                       1/6  , 1/6  , 0.5  , 1/6;
                       1/6  , 1/6  , 1/6  , 0.5];
            weights = [-4/5 , 9/20 , 9/20 , 9/20 , 9/20];
        otherwise
            error('No such Nq at this moment\n Choose between 1, 4, 5')            
    end
end