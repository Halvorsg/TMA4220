function [u_sol, u_exact] = oppgave2(N)
% a)
addpath ..\Grids
addpath ..\Oppgave1
theta = 3*pi/2;
% for i = 1:3
% Case 1
%     [p,tri,edge] = getSlice(N,3*pi/2);
%     figure
%     TR = triangulation(tri, p);
%     triplot(TR)
% % %     s = sprintf('N = %i triangles', N);
% % %     title(s)
% % % end
[p,tri,edge] = getSlice(N,theta);
TR = triangulation(tri, p);
f = @(x,y) 16*pi^2*x.*y.*(x.^2+y.^2).*sin(2*pi*(x.^2+y.^2))-24*x.*y.*pi.*cos(2*pi*(x.^2+y.^2));
exact_solution = @(x,y) x.*y.*sin(2*pi.*(x.^2+y.^2));
g = @(p) p(:,1).*p(:,2).*sin(2*pi.*(p(:,1).^2+p(:,2).^2));
% % Case 2
% [p,tri,edge] = getDisk(N);
% TR = triangulation(tri, p);
% g = @(x,y) 1;
% exact_solution = @(x,y) (1-x.^2-y.^2)/4;

A = spalloc(length(p),length(p),10*length(p)); F = zeros(length(p),1);
phi = [eye(2),[-1;-1]];


for i = 1:length(tri)
    p1 = p(tri(i,1),:)'; p2 = p(tri(i,2),:)'; p3 = p(tri(i,3),:)'; 
%     x = [p1,p2,p3];
    K = [p1-p3 , p2-p3];
    jacDet = abs(det(K));
%Stiffness element matrix - Ak
    G = K'\phi;
    Ak = jacDet*(G'*G)/2;
%Load elemen vector - Fk
    Fk = quadrature2D(p1,p2,p3,3,f);
    %Adding into Stiffnessmatrix and load vector
    A(tri(i,:),tri(i,:)) = A(tri(i,:),tri(i,:)) + Ak;
    F(tri(i,:)) = F(tri(i,:)) + Fk;

end
addpath C:\Users\halvo\Documents\MATLAB\TMA4220
% [RHS,LHS,RHS_star,LHS_star] = run_2d( N,g );
u_sol = zeros(size(p(:,1)));

[A,F,inner_vertices] = remove_boundary(A,F,edge);

% tic
% y = ismember(tri,edge);
% z = sum(y,2);
% x = tri(z<=1,:);
% m = ~ismember(x,edge);
% n = x(m);
% inner_vertices_2 = unique(n);
% toc

u = A\F;
u_sol(inner_vertices) = u;
u_exact = exact_solution(p(:,1),p(:,2));
figure
trimesh(tri,p(:,1),p(:,2),u_sol)
s = sprintf('Approximate solution with N = %i', N);
title(s)
figure
trimesh(tri,p(:,1),p(:,2),u_exact)
title('Exact solution')

end

