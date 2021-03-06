function [ Stiffness_matrix,Mass_matrix ] = eigen_value_hook( something)
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

                A(i,j)=-abs(det_J_t)*abs(det_J_t)*dot(shear_vector(:,i),C*shear_vector(:,j));
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
                Mass_mat((i-1)*3+1,(j-1)*3+1)=gauss_quad_3(x_1,x_2,x_3,x_4,5,@(x,y,z) RHS_pre_phi_1(x,y,z)*RHS_pre_phi_2(x,y,z))/(-abs(det_J_t));
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

