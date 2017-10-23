function [p tet edge] = getSphere(N)
% function [p tet edge] = getSphere(N),
% 
% description:
%      generate a mesh triangulation of the unit ball
%
% arguments:
%   - N    the number of nodes in the mesh
% returns:
%   - p     nodal points. (x,y,z)-coordinates for point i given in row i.
%   - tet   elements. Index to the four corners of element i given in row i.
%   - edge  index list of all triangles on the ball boundary

% author: Kjetil A. Johannessen
% last edit: September 2016

% approximating design
M = floor((sqrt(pi)*sqrt(12*N+pi)-3*pi)/4/pi);    % number of circles outward (not counting origin)
r = 0:1/M:1;              % radius of the different circles
totArea = sum(4*pi*r.*r); % total amount of area to distribute points on
alpha   = floor(4*pi*r.*r/totArea*N); % number of DOF in each circle

% fine tuning to get the right amount of DOF
alpha(1) = 1;

i = 2;
while(sum(alpha)>N),
	if(alpha(i)>0),
		alpha(i) = alpha(i) - 1;
	end
	i = i+1;
	if(sum(alpha(2:M)))==0,
		i = M+1;
	elseif(i>M),
		i = 2;
	end
end
while(sum(alpha)<N),
	alpha(end) = alpha(end) + 1;
end

p = [0,0,0];

k = 2;
for i=2:M+1,
	p = [p; shell(alpha(i), r(i), 2.0*pi/(M+1)*(i-1))];
end

tet      = delaunay(p(:,1), p(:,2), p(:,3));
topology = TriRep(tet, p);
edge     = freeBoundary(topology);
