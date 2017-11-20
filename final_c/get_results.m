function [V,D,p,real_inner_vertices] = get_results()

addpath ('/home/shomeb/h/halvorsg/Matlab/MATLAB/TMA4220/Project 2/for charles/Grids');
addpath ('/home/shomeb/h/halvorsg/Matlab/MATLAB/TMA4220/Project 2/for charles/Oppgave1');
%% Get triangulation
%[p,tri,edge] = getDefinedBox(-1,1,-1,1,-1,1,N);
%TR = triangulation(tri,p);
%tetramesh(TR)
[p,tri] = getMesh('plate_msh_coarser.msh');
%TR = triangulation(tri,p);
%tetramesh(TR)

%[p,tri] = getMesh('rectangle.msh');
C  = @(E,v)reshape([(E.*(v-1.0))./(v+v.^2.*2.0-1.0),-(E.*v)./(v+v.^2.*2.0-1.0),-(E.*v)./(v+v.^2.*2.0-1.0),0.0,0.0,0.0,-(E.*v)./(v+v.^2.*2.0-1.0),(E.*(v-1.0))./(v+v.^2.*2.0-1.0),-(E.*v)./(v+v.^2.*2.0-1.0),0.0,0.0,0.0,-(E.*v)./(v+v.^2.*2.0-1.0),-(E.*v)./(v+v.^2.*2.0-1.0),(E.*(v-1.0))./(v+v.^2.*2.0-1.0),0.0,0.0,0.0,0.0,0.0,0.0,(E.*(1.0./2.0))./(v+1.0),0.0,0.0,0.0,0.0,0.0,0.0,(E.*(1.0./2.0))./(v+1.0),0.0,0.0,0.0,0.0,0.0,0.0,(E.*(1.0./2.0))./(v+1.0)],[6,6]);
B  = @(p1,p2,p3,p4)reshape([p2(2).*p3(3)-p3(2).*p2(3)-p2(2).*p4(3)+p4(2).*p2(3)+p3(2).*p4(3)-p4(2).*p3(3),0.0,0.0,0.0,p2(1).*p3(2)-p3(1).*p2(2)-p2(1).*p4(2)+p4(1).*p2(2)+p3(1).*p4(2)-p4(1).*p3(2),-p2(1).*p3(3)+p3(1).*p2(3)+p2(1).*p4(3)-p4(1).*p2(3)-p3(1).*p4(3)+p4(1).*p3(3),0.0,-p2(1).*p3(3)+p3(1).*p2(3)+p2(1).*p4(3)-p4(1).*p2(3)-p3(1).*p4(3)+p4(1).*p3(3),0.0,p2(1).*p3(2)-p3(1).*p2(2)-p2(1).*p4(2)+p4(1).*p2(2)+p3(1).*p4(2)-p4(1).*p3(2),0.0,p2(2).*p3(3)-p3(2).*p2(3)-p2(2).*p4(3)+p4(2).*p2(3)+p3(2).*p4(3)-p4(2).*p3(3),0.0,0.0,p2(1).*p3(2)-p3(1).*p2(2)-p2(1).*p4(2)+p4(1).*p2(2)+p3(1).*p4(2)-p4(1).*p3(2),-p2(1).*p3(3)+p3(1).*p2(3)+p2(1).*p4(3)-p4(1).*p2(3)-p3(1).*p4(3)+p4(1).*p3(3),p2(2).*p3(3)-p3(2).*p2(3)-p2(2).*p4(3)+p4(2).*p2(3)+p3(2).*p4(3)-p4(2).*p3(3),0.0,-p1(2).*p3(3)+p3(2).*p1(3)+p1(2).*p4(3)-p4(2).*p1(3)-p3(2).*p4(3)+p4(2).*p3(3),0.0,0.0,0.0,-p1(1).*p3(2)+p3(1).*p1(2)+p1(1).*p4(2)-p4(1).*p1(2)-p3(1).*p4(2)+p4(1).*p3(2),p1(1).*p3(3)-p3(1).*p1(3)-p1(1).*p4(3)+p4(1).*p1(3)+p3(1).*p4(3)-p4(1).*p3(3),0.0,p1(1).*p3(3)-p3(1).*p1(3)-p1(1).*p4(3)+p4(1).*p1(3)+p3(1).*p4(3)-p4(1).*p3(3),0.0,-p1(1).*p3(2)+p3(1).*p1(2)+p1(1).*p4(2)-p4(1).*p1(2)-p3(1).*p4(2)+p4(1).*p3(2),0.0,-p1(2).*p3(3)+p3(2).*p1(3)+p1(2).*p4(3)-p4(2).*p1(3)-p3(2).*p4(3)+p4(2).*p3(3),0.0,0.0,-p1(1).*p3(2)+p3(1).*p1(2)+p1(1).*p4(2)-p4(1).*p1(2)-p3(1).*p4(2)+p4(1).*p3(2),p1(1).*p3(3)-p3(1).*p1(3)-p1(1).*p4(3)+p4(1).*p1(3)+p3(1).*p4(3)-p4(1).*p3(3),-p1(2).*p3(3)+p3(2).*p1(3)+p1(2).*p4(3)-p4(2).*p1(3)-p3(2).*p4(3)+p4(2).*p3(3),0.0,p1(2).*p2(3)-p2(2).*p1(3)-p1(2).*p4(3)+p4(2).*p1(3)+p2(2).*p4(3)-p4(2).*p2(3),0.0,0.0,0.0,p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p4(2)+p4(1).*p1(2)+p2(1).*p4(2)-p4(1).*p2(2),-p1(1).*p2(3)+p2(1).*p1(3)+p1(1).*p4(3)-p4(1).*p1(3)-p2(1).*p4(3)+p4(1).*p2(3),0.0,-p1(1).*p2(3)+p2(1).*p1(3)+p1(1).*p4(3)-p4(1).*p1(3)-p2(1).*p4(3)+p4(1).*p2(3),0.0,p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p4(2)+p4(1).*p1(2)+p2(1).*p4(2)-p4(1).*p2(2),0.0,p1(2).*p2(3)-p2(2).*p1(3)-p1(2).*p4(3)+p4(2).*p1(3)+p2(2).*p4(3)-p4(2).*p2(3),0.0,0.0,p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p4(2)+p4(1).*p1(2)+p2(1).*p4(2)-p4(1).*p2(2),-p1(1).*p2(3)+p2(1).*p1(3)+p1(1).*p4(3)-p4(1).*p1(3)-p2(1).*p4(3)+p4(1).*p2(3),p1(2).*p2(3)-p2(2).*p1(3)-p1(2).*p4(3)+p4(2).*p1(3)+p2(2).*p4(3)-p4(2).*p2(3),0.0,-p1(2).*p2(3)+p2(2).*p1(3)+p1(2).*p3(3)-p3(2).*p1(3)-p2(2).*p3(3)+p3(2).*p2(3),0.0,0.0,0.0,-p1(1).*p2(2)+p2(1).*p1(2)+p1(1).*p3(2)-p3(1).*p1(2)-p2(1).*p3(2)+p3(1).*p2(2),p1(1).*p2(3)-p2(1).*p1(3)-p1(1).*p3(3)+p3(1).*p1(3)+p2(1).*p3(3)-p3(1).*p2(3),0.0,p1(1).*p2(3)-p2(1).*p1(3)-p1(1).*p3(3)+p3(1).*p1(3)+p2(1).*p3(3)-p3(1).*p2(3),0.0,-p1(1).*p2(2)+p2(1).*p1(2)+p1(1).*p3(2)-p3(1).*p1(2)-p2(1).*p3(2)+p3(1).*p2(2),0.0,-p1(2).*p2(3)+p2(2).*p1(3)+p1(2).*p3(3)-p3(2).*p1(3)-p2(2).*p3(3)+p3(2).*p2(3),0.0,0.0,-p1(1).*p2(2)+p2(1).*p1(2)+p1(1).*p3(2)-p3(1).*p1(2)-p2(1).*p3(2)+p3(1).*p2(2),p1(1).*p2(3)-p2(1).*p1(3)-p1(1).*p3(3)+p3(1).*p1(3)+p2(1).*p3(3)-p3(1).*p2(3),-p1(2).*p2(3)+p2(2).*p1(3)+p1(2).*p3(3)-p3(2).*p1(3)-p2(2).*p3(3)+p3(2).*p2(3),0.0],[6,12]);
A6 = @(p1,p2,p3,p4)p1(1).*p2(2).*p3(3)-p1(1).*p3(2).*p2(3)-p2(1).*p1(2).*p3(3)+p2(1).*p3(2).*p1(3)+p3(1).*p1(2).*p2(3)-p3(1).*p2(2).*p1(3)-p1(1).*p2(2).*p4(3)+p1(1).*p4(2).*p2(3)+p2(1).*p1(2).*p4(3)-p2(1).*p4(2).*p1(3)-p4(1).*p1(2).*p2(3)+p4(1).*p2(2).*p1(3)+p1(1).*p3(2).*p4(3)-p1(1).*p4(2).*p3(3)-p3(1).*p1(2).*p4(3)+p3(1).*p4(2).*p1(3)+p4(1).*p1(2).*p3(3)-p4(1).*p3(2).*p1(3)-p2(1).*p3(2).*p4(3)+p2(1).*p4(2).*p3(3)+p3(1).*p2(2).*p4(3)-p3(1).*p4(2).*p2(3)-p4(1).*p2(2).*p3(3)+p4(1).*p3(2).*p2(3);

rho = 1;
e = rho*ones(12,1);
Mk = @(p1,p2,p3,p4) 1/(120*A6(p1,p2,p3,p4))*spdiags([e,e,e,2.*e,e,e,e],-9:3:9,12,12);
%A6 is determinant <==> Volume*6 for tetrahedra
fcn = @plus; %For parallellization




E = 235000;
v = 0.346;
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
parfor i = 1:length(tri)
    i
    %     p1 = p(tri(i,1),:)'; p2 = p(tri(i,2),:)'; p3 = p(tri(i,3),:)'; p4 = p(tri(i,4),:)'; 
    %Finding the right indices for A by i = 3*(�-1) + d 
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
%y = (~ismember(tri,edge)).*tri;
%inner_vertices_2 = unique(y);
%inner_vertices_2 = inner_vertices_2(2:end);

%tempIV = [3*(inner_vertices_2-1)+1 ; 3*(inner_vertices_2-1)+2 ; 3*(inner_vertices_2-1)+3]';
%inner_vertices = tempIV(:);
inner_vertices=zeros(length(p)*3,1);
boundary_points=zeros(length(p)*3,1);
step=1;
b_step=1;
acc=0.001;
for iterator=1:length(p)
%if p(1,1)~=p(iterator,1) || p(1,2)~=p(iterator,2) || p(1,3)~=p(iterator,3)
if (p(iterator,1)>(0+acc) && p(iterator,1)<(100-acc) && p(iterator,2)>(0+acc) && p(iterator,2)<(100-acc)) %|| ((p(iterator,3)>0.27+acc || p(iterator,3)<0.23)) && (p(iterator,1)<(0+acc) || (p(iterator,1)>(10-acc) || p(iterator,2)<(0+acc) || p(iterator,2)>(10-acc)))
    inner_vertices(step)=iterator;
    step=step+1;
else
    boundary_points(b_step)=iterator;
    b_step=b_step+1;
end
end

real_inner_vertices=inner_vertices(1:step-1);
boundary_points=boundary_points(1:b_step-1);



%% xyz lock
%{

step_h=3*(step-1);
test_p=zeros(length(p),3);
test_p(real_inner_vertices,:)=p(real_inner_vertices,:);

inner_vertices(1:3:step_h)=real_inner_vertices*3-2;
inner_vertices(2:3:step_h)=real_inner_vertices*3-1;
inner_vertices(3:3:step_h)=real_inner_vertices*3;
inner_vertices=inner_vertices(1:step_h);
%}

%% z lock only
%{
step_h=3*(step-1);

inner_vertices(1:3:step_h)=real_inner_vertices*3-2;
inner_vertices(2:3:step_h)=real_inner_vertices*3-1;
inner_vertices(3:3:step_h)=real_inner_vertices*3;
%inner_vertices=inner_vertices(1:step_h);
xy_b_p_h=length(p)-(step-1);
xy_b_p_h=xy_b_p_h*2+step_h;
inner_vertices(step_h+1:2:xy_b_p_h)=boundary_points*3-2;
inner_vertices(step_h+2:2:xy_b_p_h)=boundary_points*3-1;
inner_vertices=inner_vertices(1:xy_b_p_h);
%}

%% xy_lock_only

step_h=3*(step-1);

inner_vertices(1:3:step_h)=real_inner_vertices*3-2;
inner_vertices(2:3:step_h)=real_inner_vertices*3-1;
inner_vertices(3:3:step_h)=real_inner_vertices*3;
%inner_vertices=inner_vertices(1:step_h);
xy_b_p_h=length(p)-(step-1);
xy_b_p_h=xy_b_p_h+step_h;
inner_vertices(step_h+1:xy_b_p_h)=boundary_points*3;
inner_vertices=inner_vertices(1:xy_b_p_h);

 A = A(inner_vertices,inner_vertices);
 M = M(inner_vertices,inner_vertices);
%F = F(inner_vertices);
A = sparse(A);

%% We now have M and A, want to find the eigenvalues inv(M)*A
tic
[V,D] = eigs(M,A,50);
[D,index] = sort(diag(D)); %Eigenvalues sorted from largest to smallest
V = V(:,index); % Eigenvectors corresponding to eigenvalues in D
eigen = toc;
fprintf('Eigenvalues = %f\n',eigen)

%% surface contour mode


%{
time=(2/maxiter)*j*pi;
        for i=1:length(real_inner_vertices)
        pre_plot_vec(i,:)=p(real_inner_vertices(i),:)+scaling*[V(i*3-2,1),V(i*3-1,1),V(i*3,1)]*sin(time);
        end
        %for i=1:length(boundary_points)
        %pre_plot_vec(length(real_inner_vertices)+i,:)=p(boundary_points(i),:)+scaling*[V(length(real_inner_vertices)*3+i*2-1,1),V(length(real_inner_vertices)*3+i*2,1),0]*sin(time);
        %end
  %}      
scaling = 10;

x=zeros(length(p),1);
y=zeros(length(p),1);
z=zeros(length(p),1);

step=1;
  for i=1:length(real_inner_vertices)
    if p(real_inner_vertices(i),3)>0.48
    x(step)=p(real_inner_vertices(i),1)+V(i*3-2,12)*scaling;
    y(step)=p(real_inner_vertices(i),2)+V(i*3-1,12)*scaling;
    z(step)=p(real_inner_vertices(i),3)+V(i*3,12)*scaling;
    step=step+1;
    end
  end
    x=x(1:step-1);
    y=y(1:step-1);
    z=z(1:step-1);
   
%tri_2d=delaunay(x,y);
%trisurf(tri_2d,x,y,z);
[m_x,m_y] = meshgrid(0:1:100, 0:1:100);
test_plot = griddata(x,y,z,m_x,m_y);%   
contour(m_x,m_y,test_plot);
Animation=figure(1);

%% surface animation mode
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

%% tetramesh animation mode
%{
scaling = 7;
f = figure('visible', 'off');
%pre_plot_vec = zeros(length(real_inner_vertices)+length(boundary_points),3);
pre_plot_vec = zeros(length(real_inner_vertices),3);

maxiter = 40;

for j=1:maxiter
    j  
    step=1;
        time=(2/maxiter)*j*pi;
        for i=1:length(real_inner_vertices)
        pre_plot_vec(i,:)=p(real_inner_vertices(i),:)+scaling*[V(i*3-2,6),V(i*3-1,6),V(i*3,6)]*sin(time);
        end
        %for i=1:length(boundary_points)
        %pre_plot_vec(length(real_inner_vertices)+i,:)=p(boundary_points(i),:)+scaling*[V(length(real_inner_vertices)*3+i*2-1,1),V(length(real_inner_vertices)*3+i*2,1),0]*sin(time);
        %end
        
        plot_vec=p;
        %plot_vec([real_inner_vertices;boundary_points],:)=pre_plot_vec(:,:);
          plot_vec(real_inner_vertices,:)=pre_plot_vec(:,:);

    f=figure('visible', 'off');
    TR = triangulation(tri,plot_vec);
    tetramesh(TR);

    Animation(j)=getframe(f);
end

f=figure('visible','on')
%movie(gcf,Animation,20)
%}

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





