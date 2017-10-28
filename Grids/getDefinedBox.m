function [p tet edge] = getDefinedBox(x1,x2,y1,y2,z1,z2, N)
% function [p tet edge] = getBox(N),
% 
% description:
%      generate a mesh triangulation of the unit box (-1,1)^3
%
% arguments:
%   - N    the number of nodes in each spatial direction (N^3 total nodes)
% returns:
%   - p     nodal points. (x,y,z)-coordinates for point i given in row i.
%   - tet   elements. Index to the four corners of element i given in row i.
%   - edge  index list of all nodal points on the outer edge (r=1)

% author: Kjetil A. Johannessen
% last edit: September 2016
X = max([abs(x1-x2),abs(y1-y2),abs(z1-z2)]);
N1 = max(ceil(N*abs(x1-x2)/X),3);
N2 = max(ceil(N*abs(y1-y2)/X),3);
N3 = max(ceil(N*abs(z1-z2)/X),3);

[y x z] = meshgrid(linspace(x1,x2,N1), linspace(y1,y2,N2), linspace(z1,z2,N3));
p = [x(:), y(:), z(:)];

tet  = delaunay(p(:,1), p(:,2), p(:,3));
topology = TriRep(tet, p);
edge     = freeBoundary(topology);
