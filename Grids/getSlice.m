function [p,tri,edge] = getSlice(N, theta)
% function [p tri edge] = getSlice(N, theta),
% 
% description:
%      generate a mesh triangulation of the unit disk
%
% arguments:
%   - N     the number of nodes in the mesh
%   - theta angle of slice
% returns:
%   - p      nodal points. (x,y)-coordinates for point i given in row i.
%   - tri    elements. Index to the three corners of element i given in row i.
%   - edge   edge lines. Index list to the two corners of edge line i given in row i

% author: Kjetil A. Johannessen
% last edit: September 2016

% approximating design
M = floor(sqrt(2*N/theta));    % number of circles outward (not counting origin)
r = 0:1/M:1;              % radius of the different circles
alpha = floor(theta*r*M);  % number of DOF in each circle

if N<3
  throw(MException('getSlice', 'Error: Too small input. Need meshsize bigger than 2'));
end
  

% fine tuning to get the right amount of DOF
alpha(1) = 1;

i = M+1;
while(sum(alpha)>N),
	if(alpha(i)>0),
		alpha(i) = alpha(i) - 1;
	end
	i = i-1;
	if(i<2)
		i = M+1;
	end
end
while(sum(alpha)<N),
	alpha(end) = alpha(end) + 1;
end

% special case small examples
if M == 1 || (M == 2 && sum(alpha(2))<3)
  i = 0:N-2;
  t = theta/(N-2) * i';
  p = [0,0;cos(t), sin(t)];
  edge = [2:N-1;3:N];
  tri  = [ones(N-2,1), edge'];
  edge = [edge'; N,1; 1,2];
  return;
end

p = zeros(N,2);
p(1,:) = [0,0];
edge = zeros(numel(alpha)*2-2,2);

k = 2;
for i=2:M+1,
	t = 0;
	for j=1:alpha(i),
		p(k,:) = [cos(t), sin(t)]*r(i);
		t = t + theta/(alpha(i)-1);
		k = k+1;
	end
  edge(2*(i-2)+1:2*(i-1),:) = [k-1, k-1-alpha(i); k-alpha(i)-alpha(i-1), k-alpha(i)];
end
k = N-alpha(end)+1;
edge = [edge; [k:N-1]', [k+1:N]'];

force = [N, k-2; k, k-alpha(end-1)+1];
tri   = delaunayTriangulation(p(:,1), p(:,2), force);
tri   = tri.ConnectivityList;

if theta > pi % cut out non-convex part
  e = unique(edge(:));
  j = zeros((numel(alpha)-1)*2-1,1);
  k = 1;
  for i=1:size(tri,1)
    if numel(intersect(tri(i,:), e)) == 3
      j(k) = i;
      k    = k + 1;
    end
  end
  tri(j,:) = [];
end
