function oneDimensional

addpath C:\Users\halvo\Documents\MATLAB\TMA4220\Project1\Grids
addpath C:\Users\halvo\Documents\MATLAB\TMA4220\Project1\Oppgave1
addpath C:\Users\halvo\Documents\MATLAB\Elementmetoden

N = 100;
x = linspace(0,pi);
g = @(x) sin(x);

exact_solution = @(x,y) -sin(x);

A = spalloc(N,N,N*3); F = zeros(N,1);
phi = [1,-1];
H = @(x,y) g(x,y)*[x;1-x];

for i = 1:length(N)
    p1 = x(i); p2 = x(i+1); 
    J = [p1-p2];
    jacDet = abs(det(J));
%Stiffness element matrix - Ak
    G = J\phi;
    Ak = jacDet*(G'*G);
%Load elemen vector - Fk
%     I = element_quadrature(Nq,g,xk,hk)
%     Fk = quadrature2D(p1,p2,p3,4,H);
% %Adding into Stiffnessmatrix and load vector
%     A(tri(i,:),tri(i,:)) = A(tri(i,:),tri(i,:)) + Ak;
%     F(tri(i,:)) = F(tri(i,:)) + Fk;

end
u_sol = zeros(size(p(:,1)));
[A,F,inner_vertices] = remove_boundary(A,F,edge);

u = A\F;
u_sol(inner_vertices) = u;
u = exact_solution(p(:,1),p(:,2));
figure
trimesh(tri,p(:,1),p(:,2),u_sol)
figure
trimesh(tri,p(:,1),p(:,2),u)
errors(j) = max(abs(u_sol - u));
end


