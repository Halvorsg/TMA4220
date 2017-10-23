function [p tri edge] = getDisc(N),
% function [p tri edge] = getDisc(N),
% 
% description:
%      generate a mesh triangulation of the unit disk
%
% arguments:
%   - N    the number of nodes in the mesh
% returns:
%   - p     nodal points. (x,y)-coordinates for point i given in row i.
%   - tri   elements. Index to the three corners of element i given in row i.
%   - edge  edge lines. Index list to the two corners of edge line i given in row i

% author: Kjetil A. Johannessen
% last edit: September 2012

% approximating design
M = floor(sqrt(N/pi));    % number of circles outward (not counting origin)
r = 0:1/M:1;              % radius of the different circles
alpha = floor(2*pi*r*M);  % number of DOF in each circle

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

% theta = start-angle for circle i
theta = pi./alpha;
theta(1:2:end) = 0;

p = zeros(N,2);
p(1,:) = [0,0];

k = 2;
for i=2:M+1,
	t = theta(i);
	for j=1:alpha(i),
		p(k,:) = [cos(t), sin(t)]*r(i);
		t = t + 2*pi/alpha(i);
		k = k+1;
	end
end

tri  = delaunay(p(:,1), p(:,2));
edge = N-alpha(end)+1:N;
edge = [edge', edge'+1];
edge(end) = N-alpha(end)+1;

