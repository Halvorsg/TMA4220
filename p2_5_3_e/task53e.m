function [ux,uy,uz,u_exact] = task53e(N)

addpath ..\Grids
addpath ..\Oppgave1
%% Get triangulation
[p,tri,edge] = getDefinedBox(-1,1,-1,1,-1,1,N);
TR = triangulation(tri,p);
tetramesh(TR)
% [p,tri] = getMesh('rectangle.msh');

% [p,tri] = getMesh('rectangle.msh');
C  = @(E,v)reshape([(E.*(v-1.0))./(v+v.^2.*2.0-1.0),-(E.*v)./(v+v.^2.*2.0-1.0),-(E.*v)./(v+v.^2.*2.0-1.0),0.0,0.0,0.0,-(E.*v)./(v+v.^2.*2.0-1.0),(E.*(v-1.0))./(v+v.^2.*2.0-1.0),-(E.*v)./(v+v.^2.*2.0-1.0),0.0,0.0,0.0,-(E.*v)./(v+v.^2.*2.0-1.0),-(E.*v)./(v+v.^2.*2.0-1.0),(E.*(v-1.0))./(v+v.^2.*2.0-1.0),0.0,0.0,0.0,0.0,0.0,0.0,(E.*(1.0./2.0))./(v+1.0),0.0,0.0,0.0,0.0,0.0,0.0,(E.*(1.0./2.0))./(v+1.0),0.0,0.0,0.0,0.0,0.0,0.0,(E.*(1.0./2.0))./(v+1.0)],[6,6]);
B  = @(p1,p2,p3,p4)reshape([p2(2).*p3(3)-p3(2).*p2(3)-p2(2).*p4(3)+p4(2).*p2(3)+p3(2).*p4(3)-p4(2).*p3(3),0.0,0.0,0.0,p2(1).*p3(2)-p3(1).*p2(2)-p2(1).*p4(2)+p4(1).*p2(2)+p3(1).*p4(2)-p4(1).*p3(2),-p2(1).*p3(3)+p3(1).*p2(3)+p2(1).*p4(3)-p4(1).*p2(3)-p3(1).*p4(3)+p4(1).*p3(3),0.0,-p2(1).*p3(3)+p3(1).*p2(3)+p2(1).*p4(3)-p4(1).*p2(3)-p3(1).*p4(3)+p4(1).*p3(3),0.0,p2(1).*p3(2)-p3(1).*p2(2)-p2(1).*p4(2)+p4(1).*p2(2)+p3(1).*p4(2)-p4(1).*p3(2),0.0,p2(2).*p3(3)-p3(2).*p2(3)-p2(2).*p4(3)+p4(2).*p2(3)+p3(2).*p4(3)-p4(2).*p3(3),0.0,0.0,p2(1).*p3(2)-p3(1).*p2(2)-p2(1).*p4(2)+p4(1).*p2(2)+p3(1).*p4(2)-p4(1).*p3(2),-p2(1).*p3(3)+p3(1).*p2(3)+p2(1).*p4(3)-p4(1).*p2(3)-p3(1).*p4(3)+p4(1).*p3(3),p2(2).*p3(3)-p3(2).*p2(3)-p2(2).*p4(3)+p4(2).*p2(3)+p3(2).*p4(3)-p4(2).*p3(3),0.0,-p1(2).*p3(3)+p3(2).*p1(3)+p1(2).*p4(3)-p4(2).*p1(3)-p3(2).*p4(3)+p4(2).*p3(3),0.0,0.0,0.0,-p1(1).*p3(2)+p3(1).*p1(2)+p1(1).*p4(2)-p4(1).*p1(2)-p3(1).*p4(2)+p4(1).*p3(2),p1(1).*p3(3)-p3(1).*p1(3)-p1(1).*p4(3)+p4(1).*p1(3)+p3(1).*p4(3)-p4(1).*p3(3),0.0,p1(1).*p3(3)-p3(1).*p1(3)-p1(1).*p4(3)+p4(1).*p1(3)+p3(1).*p4(3)-p4(1).*p3(3),0.0,-p1(1).*p3(2)+p3(1).*p1(2)+p1(1).*p4(2)-p4(1).*p1(2)-p3(1).*p4(2)+p4(1).*p3(2),0.0,-p1(2).*p3(3)+p3(2).*p1(3)+p1(2).*p4(3)-p4(2).*p1(3)-p3(2).*p4(3)+p4(2).*p3(3),0.0,0.0,-p1(1).*p3(2)+p3(1).*p1(2)+p1(1).*p4(2)-p4(1).*p1(2)-p3(1).*p4(2)+p4(1).*p3(2),p1(1).*p3(3)-p3(1).*p1(3)-p1(1).*p4(3)+p4(1).*p1(3)+p3(1).*p4(3)-p4(1).*p3(3),-p1(2).*p3(3)+p3(2).*p1(3)+p1(2).*p4(3)-p4(2).*p1(3)-p3(2).*p4(3)+p4(2).*p3(3),0.0,p1(2).*p2(3)-p2(2).*p1(3)-p1(2).*p4(3)+p4(2).*p1(3)+p2(2).*p4(3)-p4(2).*p2(3),0.0,0.0,0.0,p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p4(2)+p4(1).*p1(2)+p2(1).*p4(2)-p4(1).*p2(2),-p1(1).*p2(3)+p2(1).*p1(3)+p1(1).*p4(3)-p4(1).*p1(3)-p2(1).*p4(3)+p4(1).*p2(3),0.0,-p1(1).*p2(3)+p2(1).*p1(3)+p1(1).*p4(3)-p4(1).*p1(3)-p2(1).*p4(3)+p4(1).*p2(3),0.0,p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p4(2)+p4(1).*p1(2)+p2(1).*p4(2)-p4(1).*p2(2),0.0,p1(2).*p2(3)-p2(2).*p1(3)-p1(2).*p4(3)+p4(2).*p1(3)+p2(2).*p4(3)-p4(2).*p2(3),0.0,0.0,p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p4(2)+p4(1).*p1(2)+p2(1).*p4(2)-p4(1).*p2(2),-p1(1).*p2(3)+p2(1).*p1(3)+p1(1).*p4(3)-p4(1).*p1(3)-p2(1).*p4(3)+p4(1).*p2(3),p1(2).*p2(3)-p2(2).*p1(3)-p1(2).*p4(3)+p4(2).*p1(3)+p2(2).*p4(3)-p4(2).*p2(3),0.0,-p1(2).*p2(3)+p2(2).*p1(3)+p1(2).*p3(3)-p3(2).*p1(3)-p2(2).*p3(3)+p3(2).*p2(3),0.0,0.0,0.0,-p1(1).*p2(2)+p2(1).*p1(2)+p1(1).*p3(2)-p3(1).*p1(2)-p2(1).*p3(2)+p3(1).*p2(2),p1(1).*p2(3)-p2(1).*p1(3)-p1(1).*p3(3)+p3(1).*p1(3)+p2(1).*p3(3)-p3(1).*p2(3),0.0,p1(1).*p2(3)-p2(1).*p1(3)-p1(1).*p3(3)+p3(1).*p1(3)+p2(1).*p3(3)-p3(1).*p2(3),0.0,-p1(1).*p2(2)+p2(1).*p1(2)+p1(1).*p3(2)-p3(1).*p1(2)-p2(1).*p3(2)+p3(1).*p2(2),0.0,-p1(2).*p2(3)+p2(2).*p1(3)+p1(2).*p3(3)-p3(2).*p1(3)-p2(2).*p3(3)+p3(2).*p2(3),0.0,0.0,-p1(1).*p2(2)+p2(1).*p1(2)+p1(1).*p3(2)-p3(1).*p1(2)-p2(1).*p3(2)+p3(1).*p2(2),p1(1).*p2(3)-p2(1).*p1(3)-p1(1).*p3(3)+p3(1).*p1(3)+p2(1).*p3(3)-p3(1).*p2(3),-p1(2).*p2(3)+p2(2).*p1(3)+p1(2).*p3(3)-p3(2).*p1(3)-p2(2).*p3(3)+p3(2).*p2(3),0.0],[6,12]);
A6 = @(p1,p2,p3,p4)p1(1).*p2(2).*p3(3)-p1(1).*p3(2).*p2(3)-p2(1).*p1(2).*p3(3)+p2(1).*p3(2).*p1(3)+p3(1).*p1(2).*p2(3)-p3(1).*p2(2).*p1(3)-p1(1).*p2(2).*p4(3)+p1(1).*p4(2).*p2(3)+p2(1).*p1(2).*p4(3)-p2(1).*p4(2).*p1(3)-p4(1).*p1(2).*p2(3)+p4(1).*p2(2).*p1(3)+p1(1).*p3(2).*p4(3)-p1(1).*p4(2).*p3(3)-p3(1).*p1(2).*p4(3)+p3(1).*p4(2).*p1(3)+p4(1).*p1(2).*p3(3)-p4(1).*p3(2).*p1(3)-p2(1).*p3(2).*p4(3)+p2(1).*p4(2).*p3(3)+p3(1).*p2(2).*p4(3)-p3(1).*p4(2).*p2(3)-p4(1).*p2(2).*p3(3)+p4(1).*p3(2).*p2(3);

rho = 1;
e = rho*ones(12,1);
Mk = @(p1,p2,p3,p4) 1/(120*A6(p1,p2,p3,p4))*spdiags([e,e,e,2.*e,e,e,e],-9:3:9,12,12);
%A6 is determinant <==> Volume*6 for tetrahedra
fcn = @plus; %For parallellization




E = 1;
v = 0.3;
fx = @(x,y,z) (E/((1+v)*(1-2*v)))* ( 4*v*x*y*(z^2 - 1) - 2*(y^2 - 1)*(v - 1/2)*(x^2 + 2*z*x - 1) - 2*(z^2 - 1)*(v - 1/2)*(x^2 + 2*y*x - 1) - 2*(y^2 - 1)*(z^2 - 1)*(v - 1) + 4*v*x*z*(y^2 - 1));
fy = @(x,y,z) (E/((1+v)*(1-2*v)))* ( 4*v*x*y*(z^2 - 1) - 2*(x^2 - 1)*(v - 1/2)*(y^2 + 2*z*y - 1) - 2*(z^2 - 1)*(v - 1/2)*(y^2 + 2*x*y - 1) - 2*(x^2 - 1)*(z^2 - 1)*(v - 1) + 4*v*y*z*(x^2 - 1));
fz = @(x,y,z) (E/((1+v)*(1-2*v)))* ( 4*v*x*z*(y^2 - 1) - 2*(x^2 - 1)*(v - 1/2)*(z^2 + 2*y*z - 1) - 2*(y^2 - 1)*(v - 1/2)*(z^2 + 2*x*z - 1) - 2*(x^2 - 1)*(y^2 - 1)*(v - 1) + 4*v*y*z*(x^2 - 1));

exact = @(x,y,z) (x.^2-1).*(y.^2-1).*(z.^2-1);

Ak = @(p1,p2,p3,p4) B(p1,p2,p3,p4)'*C(E,v)*B(p1,p2,p3,p4)*(1/6)*1/A6(p1,p2,p3,p4);

%% Preallocate and iterate over all elements
N = max(max(tri));
A = spalloc(3*N,3*N,10*N);
M = spalloc(3*N,3*N,10*N);
p1 = p(tri(:,1),:)'; p2 = p(tri(:,2),:)'; p3 = p(tri(:,3),:)'; p4 = p(tri(:,4),:)';
tic
for i = 1:length(tri)
%     p1 = p(tri(i,1),:)'; p2 = p(tri(i,2),:)'; p3 = p(tri(i,3),:)'; p4 = p(tri(i,4),:)'; 
    %Finding the right indices for A by i = 3*(î-1) + d 
    tempIndices = [3*(tri(i,:)-1)+1 ; 3*(tri(i,:)-1)+2 ; 3*(tri(i,:)-1)+3];
    indices = tempIndices(:);
    [rowMap,colMap] = get_node_mapping_matrix(indices);
    A = fcn(A,sparse(rowMap,colMap,Ak(p1(:,i),p2(:,i),p3(:,i),p4(:,i)),3*N,3*N));
    M = fcn(M,sparse(rowMap,colMap,Mk(p1(:,i),p2(:,i),p3(:,i),p4(:,i)),3*N,3*N));
    % Finding the quadrature for fx and fy separate, easier to use existing
    % code this way

end
mat = toc;
% fprintf('Matrix = %f \n Quadrature = %f\n',mat,quad)
F = zeros(length(A),1);
tic
for i = 1:length(tri)
%     p1 = p(tri(i,1),:)'; p2 = p(tri(i,2),:)'; p3 = p(tri(i,3),:)'; p4 = p(tri(i,4),:)'; 
    tempIndices = [3*(tri(i,:)-1)+1 ; 3*(tri(i,:)-1)+2 ; 3*(tri(i,:)-1)+3];
    indices = tempIndices(:);
    Fkx = quadrature3D(p1(:,i),p2(:,i),p3(:,i),p4(:,i),4,fx);
    Fky = quadrature3D(p1(:,i),p2(:,i),p3(:,i),p4(:,i),4,fy);
    Fkz = quadrature3D(p1(:,i),p2(:,i),p3(:,i),p4(:,i),4,fz);
    Fk = [Fkx';Fky';Fkz'];
    F(indices) = F(indices) + Fk(:);

end
quad = toc;
fprintf('Matrix = %f \nQuadrature = %f\n',mat,quad)

%% Implementing boundary conditions and solving
y = (~ismember(tri,edge)).*tri;
inner_vertices_2 = unique(y);
inner_vertices_2 = inner_vertices_2(2:end);

tempIV = [3*(inner_vertices_2-1)+1 ; 3*(inner_vertices_2-1)+2 ; 3*(inner_vertices_2-1)+3]';
inner_vertices = tempIV(:);


% A = A(inner_vertices,inner_vertices);
% M = M(inner_vertices,inner_vertices);
F = F(inner_vertices);
A = sparse(A);

%% We now have M and A, want to find the eigenvalues inv(M)*A
tic

[V,D] = eigs(M,A,20);
[D,index] = sort(diag(D)); %Eigenvalues sorted from largest to smallest
V = V(:,index); % Eigenvectors corresponding to eigenvalues in D
eigen = toc;
fprintf('Eigenvalues = %f\n',eigen)
%
%% surface mpde
%{
scaling = 10;
f = figure('visible', 'off');
x=zeros(length(p),1);
y=zeros(length(p),1);
z=zeros(length(p),1);


maxiter = 10;
for jtt=1:maxiter
  jtt  
step=1;
    time=(2/maxiter)*jtt*pi;
    for itt=1:length(p)
    if p(itt,3)>0.98
    x(step)=p(itt,1)+V(itt*3-2,5)*scaling*sin(time);
    y(step)=p(itt,2)+V(itt*3-1,5)*scaling*sin(time);
    z(step)=p(itt,3)+V(itt*3,5)*scaling*sin(time);
    step=step+1;
    end
    end
    x=x(1:step-1);
    y=y(1:step-1);
    z=z(1:step-1);
   
    
     f=figure('visible', 'off');
[m_x,m_y] = meshgrid(-1:.02:1, -1:.02:1);
test_plot = griddata(x,y,z,m_x,m_y);
surf(-1:0.02:1,-1:0.02:1,test_plot);
zlim([0 2])
Animation(jtt)=getframe(f);
%   

end
f=figure('visible','on')
movie(gcf,Animation,20)

%}
%% tetramesh mode

scaling = 10;
f = figure('visible', 'off');
plot_vec = zeros(length(p),3);
maxiter = 10;

for j=1:maxiter
    j  
    step=1;
        time=(2/maxiter)*j*pi;
        for i=1:length(p)
        plot_vec(i,:)=p(i,:)+scaling*[V(i*3-2,1),V(i*3-1,1),V(i*3,1)]*sin(time);
        end

    f=figure('visible', 'off');
    TR = triangulation(tri,plot_vec);
    tetramesh(TR);

    Animation(j)=getframe(f);
end

f=figure('visible','on')
movie(gcf,Animation,20)

%% example problem
%{
%% Solving the test problem
u = A\F;
u_sol = zeros(3*length(p),1);
u_sol(inner_vertices) = u;
ux = u_sol(1:3:end);
uy = u_sol(2:3:end);
uz = u_sol(3:3:end);

u_exact = exact(p(:,1),p(:,2),p(:,3));

% plot(ux)
% hold on
% plot(uy)
% plot(uz)
% plot(u_exact)
%% Visualization
%}

end





