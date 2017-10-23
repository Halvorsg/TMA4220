function pts = shell(N, r, alpha),
% function pts = shell(N, r),
% 
% description:
%      generate a point cloud on a sphere shell of radius r
%
% arguments:
%   - N     the number of nodes in the mesh
%   - r     the radius of the sphere shell
%   - alpha angle start offset (optional)
% returns:
%   - pts  nodal points. (x,y,z)-coordinates for point i given in row i.

% author: Kjetil A. Johannessen
% last edit: September 2012

if nargin < 3,
	alpha = 0.0;
end

inc = pi * (3.0 - sqrt(5));
off = 2.0 / N;
pts = zeros(N, 3);
for k = 0:N-1, 
	y = k * off - 1 + (off / 2);
	R = sqrt(1 - y*y);
	phi = k * inc + alpha;
	pts(k+1,:) = r*[cos(phi)*R, y, sin(phi)*R];
end