function [u_sol,u_exact] = oppgave4b1(N)
addpath C:\Users\halvo\Documents\MATLAB\TMA4220\Project1\Grids
addpath C:\Users\halvo\Documents\MATLAB\TMA4220\Project1\Oppgave1
%% Making tetrahedras
[p,tri,edge] = getSphere(N);
TR = triangulation(tri, p);
f = @(x,y,z) 12*pi^2*sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);
exact_solution = @(x,y,z) sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);
g = @(p) sin(2*pi*p(:,1)).*sin(2*pi*p(:,2)).*sin(2*pi*p(:,3));

%% Finding inner vertices
y = ismember(tri,edge);
z = sum(y,2);
x = tri(z<=2,:);
m = ~ismember(x,edge);
o = ismember(x,edge);
n = x(m);
inner_vertices_2 = unique(n);
%% Stiffness matrix and load vector
A = spalloc(length(p),length(p),10*length(p)); 
F = zeros(length(p),1);
phi = [eye(3),[-1;-1;-1]];
for i = 1:length(tri)
    p1 = p(tri(i,1),:)'; p2 = p(tri(i,2),:)'; p3 = p(tri(i,3),:)'; p4 = p(tri(i,4),:)';
    K = [p1-p4 , p2-p4, p3-p4];
    jacDet = abs(det(K));
%Stiffness element matrix - Ak
    G = K'\phi;
    Ak = jacDet*(G'*G)/6;
%Load elemen vector - Fk
    Fk = quadrature3D(p1,p2,p3,p4,4,f);
%Adding into Stiffnessmatrix and load vector
    A(tri(i,:),tri(i,:)) = A(tri(i,:),tri(i,:)) + Ak;
    F(tri(i,:)) = F(tri(i,:)) + Fk;
end

%% Lifting vector
gr = zeros(length(p),1);
edge_nodes = unique(edge);
gr(edge_nodes) = g(p(edge_nodes,:));
u_zero = zeros(size(p(:,1)));

G = F - A*gr;
A = A(inner_vertices_2,inner_vertices_2);
F = G(inner_vertices_2);


u = A\F;
u_zero(inner_vertices_2) = u;
u_sol = u_zero+gr;
u_exact = exact_solution(p(:,1),p(:,2),p(:,3));
% plot(u_exact)
% hold on
% plot(u_exact-u_sol)
end



