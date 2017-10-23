function I = quadratureLine2D_Show(a,b,Nq,g)
% Pre-conditions 
    % a = [x0,y0] Startpoint
    % b = [x1,y1] Endpoint
    % Nq = Number of points evaluated 
    % g = g(x,y) function evaluated 
% Post-conditions
    % I = lineintegral of g(x,y) from a to b
phis = @(t) [(1+t)/2 ; (1-t)/2];
dgamma = 1/2*sqrt((b(1)-a(1))^2 + (b(2)-a(2))^2);
f = @(t) g(1/2*((1-t)*b(1) + a(1)*(1+t)),1/2*((1-t)*b(2) + a(2)*(t+1)))*dgamma*phis(t);
I = 0;
[lamdas,weights] = gaussWeights(Nq);
    for i = 1:Nq
        I = I + f(lamdas(i))*weights(i);
    end

end
function [z,r]= gaussWeights(Nq)
    switch Nq
        case 1
            r = 2;
            z = 0;
        case 2
            r = [1,1]';
            z = [-sqrt(1/3),sqrt(1/3)]';
        case 3
            r = [5/9,8/9,5/9]';
            z = [-sqrt(3/5),0,sqrt(3/5)]';
        case 4
            r = 1/36*[18-sqrt(30),18+sqrt(30),18+sqrt(30)...
                18-sqrt(30)]';
            z = [3+2*sqrt(6/5),3-2*sqrt(6/5),3-2*sqrt(6/5),3+2*sqrt(6/5)]';
            z = (z./7).^(1/2);
            z = z.*[-1,-1,1,1]';
    end
end

