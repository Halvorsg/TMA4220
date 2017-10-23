function z = int_f(N,x1,x2,x3,y1,y2,y3)
% This function evaluates \iint_K f(x,y) dxdy using
% the Gaussian quadrature of order N where K is a
% triangle with vertices (x1,y1), (x2,y2) and (x3,y3).
xw = TriGaussPoints(N); % get quadrature points and weights
% calculate the area of the triangle
A=abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2.0
% find number of Gauss points
NP=length(xw(:,1));
z = 0.0;
for j = 1:NP
x = x1*(1-xw(j,1)-xw(j,2))+x2*xw(j,1)+x3*xw(j,2)
y = y1*(1-xw(j,1)-xw(j,2))+y2*xw(j,1)+y3*xw(j,2)
z = z + f(x,y)*xw(j,3);
end
z = A*z;
return
end

function z=f(x,y)
z=sin(x+y)*x;
return;
end

function xw = TriGaussPoints(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function TriGaussPoints provides the Gaussian points and weights %
% for the Gaussian quadrature of order n for the standard triangles. %
% %
% Input: n - the order of the Gaussian quadrature (n<=12) %
% %
% Output: xw - a n by 3 matrix: %
% 1st column gives the x-coordinates of points %
% 2nd column gives the y-coordinates of points %
% 3rd column gives the weights %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xw = zeros(n,3);
if (n == 1)
xw=[0.33333333333333 0.33333333333333 1.00000000000000];
elseif (n == 2)
xw=[0.16666666666667 0.16666666666667 0.33333333333333
0.16666666666667 0.66666666666667 0.33333333333333
0.66666666666667 0.16666666666667 0.33333333333333];
elseif (n == 3)
xw=[0.33333333333333 0.33333333333333 -0.56250000000000
0.20000000000000 0.20000000000000 0.52083333333333
0.20000000000000 0.60000000000000 0.52083333333333
0.60000000000000 0.20000000000000 0.52083333333333];
end
end