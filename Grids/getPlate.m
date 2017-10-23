function [p tri edge] = getPlate(N),
% function [p tri edge] = getPlate(N),
% 
% description:
%      generate a mesh triangulation of the reference square (-1,1)^2
%
% arguments:
%   - N    the number of nodes in each spatial direction (N^2 total nodes)
% returns:
%   - p     nodal points. (x,y)-coordinates for point i given in row i.
%   - tri   elements. Index to the three corners of element i given in row i.
%   - edge  index list of all nodal points on the outer edge (r=1)

% author: Kjetil A. Johannessen
% last edit: September 2016

p = zeros(N*N,2);
[y x] = meshgrid(linspace(-1,1,N), linspace(-1,1,N));
p = [x(:), y(:)];

tri  = delaunay(p(:,1), p(:,2));
south = [1:N-1; 2:N];
east  = [N:N:(N^2-N); (2*N):N:(N^2)];
west  = [(N^2-N+1):-N:N; (N^2-2*N+1):-N:1];
north = [N^2:-1:N^2-N+2; N^2-1:-1:N^2-N+1];
edge = [south'; east'; north'; west'];
