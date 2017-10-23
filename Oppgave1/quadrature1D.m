function I = quadrature1D(a,b,Nq,g)
% I = value of the integral
% a = startpoint
% b = endpoint
% Nq = number of integration points
% g = function pointer
f = @(x,a,b) g((b-a)/2.*x+(b+a)/2.*x./x);
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
I = dot(f(z,a,b),r)*(b-a)/2;
end