function [p tet edge] = getBox(N),
% function [p tet edge] = getBox(N),
% 
% description:
%      generate a mesh triangulation of the unit box (0,1)^3
%
% arguments:
%   - N    the number of nodes in each spatial direction (N^3 total nodes)
% returns:
%   - p     nodal points. (x,y,z)-coordinates for point i given in row i.
%   - tet   elements. Index to the four corners of element i given in row i.
%   - edge  index list of all nodal points on the outer edge (r=1)

% author: Kjetil A. Johannessen
% last edit: September 2016

[y x z] = meshgrid(linspace(0,1,N), linspace(0,1,N), linspace(0,1,N));
p = [x(:), y(:), z(:)];

tet  = delaunay(p(:,1), p(:,2), p(:,3));
topology = TriRep(tet, p);
edge     = freeBoundary(topology);
