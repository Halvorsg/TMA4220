function I = quadrature2D_Triangle(p1,p2,p3,Nq,g)

K = [p1-p3,p2-p3];
H = @(x,y) g( K(1,:)*[x;y] + p3(1) , K(2,:)*[x;y] + p3(2)).*[x;y;1-x-y];
detK = abs(det(K));
p = @(x,y) K*[x;y] + p3;
I = 0;
[points,weights] = GaussWeights2D(Nq);
for i = 1:Nq
    point = p(points(i,1),points(i,2));
    I = I + weights(i)*H(points(1),points(2));
end
I = I*detK;
end
function [lambdas, weights] = GaussWeights2D(Nq)
switch Nq
    case 1
        lambdas=[0.33333333333333 0.33333333333333];
        weights = 1;
    case 3
    lambdas=[1/6 1/6
            1/6 2/3
            2/3 1/6];
    weights = [1/3,1/3,1/3];
end
end
