function [u_sol,u_exact] = oppgave4b(N)
addpath C:\Users\halvo\Documents\MATLAB\TMA4220\Project1\Grids
addpath C:\Users\halvo\Documents\MATLAB\TMA4220\Project1\Oppgave1
%% Making tetrahedras
tic
[p,tri,edge] = getSphere(N);
TR = triangulation(tri, p);
toc
f = @(x,y,z) 12*pi^2*sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);
exact_solution = @(x,y,z) sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);
g = @(p) sin(2*pi*p(1)).*sin(2*pi*p(2)).*sin(2*pi*p(3));

A = spalloc(length(p),length(p),10*length(p)); 
Rg = zeros(length(p),1);
F = zeros(length(p),1);
phi = [eye(3),[-1;-1;-1]];

y = ismember(tri,edge);
z = sum(y,2);
edge_triangle = (z==3);
edgenum = zeros(1,1000);
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

    if edge_triangle(i) == true
%         edgenum(i) = true;
        edgenodes = double(ismember(tri(i,:),edge));
        edgenodes = edgenodes'*edgenodes;
        Fk = Fk - Ak.*edgenodes*[g(p1),g(p2),g(p3),g(p4)]';
%         B(tri(i,edgenodes), tri(i,edgenodes)) = B(tri(i,edgenodes), tri(i,edgenodes)) + Ak(edgenodes,edgenodes);
        edgenodes = double(ismember(tri(i,:),edge));
        Rg(tri(i,:)) = Rg(tri(i,:)) + quadrature3D(p1,p2,p3,p4,4,exact_solution).*edgenodes';
    end
    A(tri(i,:),tri(i,:)) = A(tri(i,:),tri(i,:)) + Ak;
    F(tri(i,:)) = F(tri(i,:)) + Fk;
end
u_sol = zeros(size(p(:,1)));

y = ismember(tri,edge);
z = sum(y,2);
x = tri(z<=2,:);
m = ~ismember(x,edge);
n = x(m);
inner_vertices_2 = unique(n);
A = A(inner_vertices_2,inner_vertices_2);
% B = B(inner_vertices_2,inner_vertices_2);
F = F(inner_vertices_2);


u = A\F;
u_sol(inner_vertices_2) = u;
u_sol = u_sol+Rg;
u_exact = exact_solution(p(:,1),p(:,2),p(:,3));
% plot(u_sol)
% hold on
% plot(u_exact)
end



