function [u_sol,u_exact] = oppgave4(N)
addpath C:\Users\halvo\Documents\MATLAB\TMA4220\Project1\Grids
addpath C:\Users\halvo\Documents\MATLAB\TMA4220\Project1\Oppgave1
%% Making tetrahedras
[p,tri,edge] = getBox(N);
TR = triangulation(tri, p);
g = @(x,y,z) 12*pi^2*sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);
exact_solution = @(x,y,z) sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);

A = spalloc(length(p),length(p),10*length(p)); 
F = zeros(length(p),1);
phi = [eye(3),[-1;-1;-1]];

tic
parfor i = 1:length(tri)
    p1 = p(tri(i,1),:)'; p2 = p(tri(i,2),:)'; p3 = p(tri(i,3),:)'; p4 = p(tri(i,4),:)';
    K = [p1-p4 , p2-p4, p3-p4];
    jacDet = abs(det(K));
%Stiffness element matrix - Ak
    G = K'\phi;
    Ak{i} = jacDet*(G'*G)/6;
%Load elemen vector - Fk
    Fk{i} = quadrature3D(p1,p2,p3,p4,4,g);
end
for i = 1:length(Fk)
%Adding into Stiffnessmatrix and load vector
    A(tri(i,:),tri(i,:)) = A(tri(i,:),tri(i,:)) + Ak{i};
    F(tri(i,:)) = F(tri(i,:)) + Fk{i};
end
toc
u_sol = zeros(size(p(:,1)));


y = ismember(tri,edge);
z = sum(y,2);
x = tri(z<=2,:);
m = ~ismember(x,edge);
n = x(m);
inner_vertices_2 = unique(n);
A = A(inner_vertices_2,inner_vertices_2);
F = F(inner_vertices_2);


u = A\F;
u_sol(inner_vertices_2) = u;
u_exact = exact_solution(p(:,1),p(:,2),p(:,3));

end



