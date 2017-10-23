function [u_sol, u_exact] = oppgave4b_2D(N)
% a)
addpath C:\Users\halvo\Documents\MATLAB\TMA4220\Project1\Grids
addpath C:\Users\halvo\Documents\MATLAB\TMA4220\Project1\Oppgave1

[p,tri,edge] = getDisk(N);
TR = triangulation(tri, p);
f = @(x,y) 8*pi^2*sin(2*pi*x).*sin(2*pi*y);
exact_solution = @(x,y) sin(2*pi*x).*sin(2*pi*y);
g = @(p) sin(2*pi*p(:,1)).*sin(2*pi*p(:,2));

A = spalloc(length(p),length(p),10*length(p)); F = zeros(length(p),1);
phi = [eye(2),[-1;-1]];

% tic
parfor i = 1:length(tri)
    p1 = p(tri(i,1),:)'; p2 = p(tri(i,2),:)'; p3 = p(tri(i,3),:)';
    K = [p1-p3 , p2-p3];
    jacDet = abs(det(K));
%Stiffness element matrix - Ak
    G = K'\phi;
    Ak{i} = jacDet*(G'*G)/2;
%Load elemen vector - Fk
    Fk{i} = quadrature2D(p1,p2,p3,4,f);
end
% Adding into Stiffnessmatrix and load vector
for i = 1:length(Fk)
    A(tri(i,:),tri(i,:)) = A(tri(i,:),tri(i,:)) + Ak{i}; 
    F(tri(i,:)) = F(tri(i,:)) + Fk{i};
end
% toc








u_sol = zeros(size(p(:,1)));


% Inner vertices
y = (~ismember(tri,edge)).*tri;
inner_vertices_2 = unique(y);
inner_vertices_2 = inner_vertices_2(2:end);

% Edge vertices
gr = zeros(length(p),1);
edge_nodes = unique(edge);
gr(edge_nodes) = g(p(edge_nodes,:));

% Lifting and solving
G = F - A*gr;
A = A(inner_vertices_2,inner_vertices_2);
F = G(inner_vertices_2);
u = A\F;

u_sol(inner_vertices_2) = u;
u_sol = u_sol+gr;
u_exact = exact_solution(p(:,1),p(:,2));
% plot(u_exact)
hold on
plot(u_sol-u_exact)
% figure
% trimesh(tri,p(:,1),p(:,2),u_sol)
% figure
% trimesh(tri,p(:,1),p(:,2),u_exact)
% figure
% trimesh(tri,p(:,1),p(:,2),u_exact-u_sol)
end