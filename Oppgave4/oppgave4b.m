function [u_sol,u_exact] = oppgave4b(N)
% addpath /home/shomeb/h/halvorsg/Documents/MATLAB~1/TMA4220/Project1/Grids
% addpath /home/shomeb/h/halvorsg/Documents/MATLAB~1/TMA4220/Project1/Oppgave1
% addpath Y:\Documents\MATLAB~1\TMA4220\Project1\Grids
% addpath Y:\Documents\MATLAB~1\TMA4220\Project1\Oppgave1
addpath ..\Grids
addpath ..\Oppgave1

%% Making tetrahedras
[p,tri,edge] = getSphere(N);
% TR = triangulation(tri, p);

f = @(x,y,z) 12*pi^2*sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);
exact_solution = @(x,y,z) sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);
g = @(p) sin(2*pi*p(:,1)).*sin(2*pi*p(:,2)).*sin(2*pi*p(:,3));
fcn = @plus;

A = spalloc(length(p),length(p),10*length(p)); 
F = zeros(length(p),1);
phi = [eye(3),[-1;-1;-1]];
parfor i = 1:length(tri)
    p1 = p(tri(i,1),:)'; p2 = p(tri(i,2),:)'; p3 = p(tri(i,3),:)'; p4 = p(tri(i,4),:)';
    K = [p1-p4 , p2-p4, p3-p4];
    jacDet = abs(det(K));
%Stiffness element matrix - Ak
    G = K'\phi;
    Ak = jacDet*(G'*G)/6;
%Load elemen vector - Fk
    Fk{i} = quadrature3D(p1,p2,p3,p4,4,f);
    [rowMap,colMap] = get_node_mapping_matrix(tri(i,:));
    A = fcn(A,sparse(rowMap,colMap,Ak,length(p),length(p)));
%Assembling
end

% Adding into load vector
for i = 1:length(Fk)
    F(tri(i,:)) = F(tri(i,:)) + Fk{i};
end
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
u_exact = exact_solution(p(:,1),p(:,2),p(:,3));
%% Visualization
F = TriScatteredInterp(p(:,1),p(:,2),p(:,3),u_sol);
[X,Y,Z] = meshgrid(-1:.1:1,-1:.1:1,-1:.1:1);
V = F(X,Y,Z);
figure
isosurface(X,Y,Z,V,0)
% tri = tri(edge_nodes,:);
% p = p(edge_nodes,:);
% u_sol = u_sol(edge_nodes);
% plot(u_exact)
% hold on
% figure
% plot(u_exact-u_sol)
end