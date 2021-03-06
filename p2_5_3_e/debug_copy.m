function [ux,uy,uz,u_exact] = debug_copy(N)

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
parfor i = 1:length(tri)
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

Wanna_know_mat_1=A;
Wanna_know_mat_2=M;

fprintf('Matrix = %f \nQuadrature = %f\n',mat,quad)

%% Implementing boundary conditions and solving
y = (~ismember(tri,edge)).*tri;
inner_vertices_2 = unique(y);
inner_vertices_2 = inner_vertices_2(2:end);

tempIV = [3*(inner_vertices_2-1)+1 ; 3*(inner_vertices_2-1)+2 ; 3*(inner_vertices_2-1)+3]';
inner_vertices = tempIV(:);

A = A(inner_vertices,inner_vertices);
M = M(inner_vertices,inner_vertices);
F = F(inner_vertices);
A = sparse(A);

%% We now have M and A, want to find the eigenvalues inv(M)*A
tic
[V,D] = eig(full(M\A));
[D,index] = sort(diag(D)); %Eigenvalues sorted from largest to smallest
V = V(:,index); % Eigenvectors corresponding to eigenvalues in D
eigen = toc;
fprintf('Eigenvalues = %f\n',eigen)

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

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
addpath ..\Grids
addpath ..\Oppgave1
%% Get triangulation
[points,Elements,edge] = getDefinedBox(-1,1,-1,1,-1,1,10);
TR = triangulation(Elements,points);
tetramesh(TR)

E = 1; v = 0.3;
  E_value= E;
        my_value= v;
        my_frac=(1+my_value)*(2);
        density=1;
        corr_vec=[0 0 0 1];
        basis_f=[eye(3),-ones(3,1)];

        
C_inv=blkdiag(-(my_value)*ones(3)+(my_value+1)*eye(3),my_frac*eye(3))/E_value;
C=inv(C_inv);

%[points,Elements] = getMesh('rectangle.msh');

Stiffness_matrix=spalloc(length(points)*3,length(points)*3,length(points)*100);
Mass_matrix=spalloc(length(points)*3,length(points)*3,length(points)*100);

cnt_dwn=length(Elements);


for iterator=1:length(Elements)
 cnt_dwn=cnt_dwn-1
      
        x_1=points(Elements(iterator,1),:);
        x_2=points(Elements(iterator,2),:);
        x_3=points(Elements(iterator,3),:);
        x_4=points(Elements(iterator,4),:);
       
        J_t=[transpose(x_1-x_4), transpose(x_2-x_4),transpose(x_3-x_4)];
        det_J_t=(det(J_t'));

  
        
        G=(transpose([J_t])\basis_f);
        small_G=G;
        temp=blkdiag(G,G,G);
        G=[temp(:,1),temp(:,5),temp(:,9),temp(:,2),temp(:,6),temp(:,10),temp(:,3),temp(:,7),temp(:,11),temp(:,4),temp(:,8),temp(:,12)];
        
        e_xx=G(1,:);
        e_yy=G(5,:);
        e_zz=G(9,:);
        
        e_xy=G(2,:)+G(4,:);
        e_xz=G(3,:)+G(7,:);
        e_yz=G(6,:)+G(8,:);
        
        shear_vector=[e_xx;e_yy;e_zz;e_xy;e_xz;e_yz];
     
     
        A=zeros(12);
      
        for i=1:12
            for j=1:12

                A(i,j)=-abs(det_J_t)*dot(shear_vector(:,i),C*shear_vector(:,j));
                %Elements(iterator,ceil(i/3))*3-mod(i+1,3)
                Stiffness_matrix(Elements(iterator,ceil(i/3))*3-mod(i+1,3),Elements(iterator,ceil(j/3))*3-mod(j+1,3))=A(i,j);

            end
        end
       
     
     
        Mass_mat=zeros(12);
        G=small_G;        
        for i=1:4
                RHS_pre_phi_2=@(x,y,z) ([x,y,z]-x_4)*G(1:3,1:3)*basis_f(:,i)+corr_vec(i);
                      
            for j=1:4
                RHS_pre_phi_1 =@(x,y,z) ([x,y,z]-x_4)*G(1:3,1:3)*basis_f(:,j)+corr_vec(j);
                Mass_mat((i-1)*3+1,(j-1)*3+1)=gauss_quad_3(x_1,x_2,x_3,x_4,5,@(x,y,z) RHS_pre_phi_1(x,y,z)*RHS_pre_phi_2(x,y,z))/(-1*det_J_t^2);
                Mass_mat((i-1)*3+2,(j-1)*3+2)=Mass_mat((i-1)*3+1,(j-1)*3+1);
                Mass_mat((i-1)*3+3,(j-1)*3+3)=Mass_mat((i-1)*3+1,(j-1)*3+1);
            end
        end
        
        for i=1:12
            for j=1:12
                Mass_matrix(Elements(iterator,ceil(i/3))*3-mod(i+1,3),Elements(iterator,ceil(j/3))*3-mod(j+1,3))=Mass_mat(i,j);
            end
        end
     
       
        
        
end
nnz(Mass_matrix)
nnz(Stiffness_matrix)
Comp_mat=Mass_matrix\Stiffness_matrix;

size(Comp_mat)



end





