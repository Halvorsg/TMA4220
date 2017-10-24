function [ soething_else ] = eigen_value_hook( something)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here




        E_value= 3.43645;
        my_value= 0.456;
        my_frac=(1-my_value)/(2)
        density=1;
        corr_vec=[0 0 0 1];
        basis_f=[eye(3),-ones(3,1)];
 %{       
   for iterator=1:length(tet)
   
        x_1=p(tet(iterator,1),:);
        x_2=p(tet(iterator,2),:);
        x_3=p(tet(iterator,3),:);
        x_4=p(tet(iterator,4),:);
      %}  
        J_t=[transpose(x_1-x_4), transpose(x_2-x_4),transpose(x_3-x_4)];
        det_J_t=det(J_t);

        basis_f=[eye(3),-ones(3,1)];
        
        G=transpose([J_t])\basis_f;
        small_G=G;
        temp=blkdiag(G,G,G);
        G=[temp(:,1),temp(:,5),temp(:,9),temp(:,2),temp(:,6),temp(:,10),temp(:,3),temp(:,7),temp(:,11),temp(:,4),temp(:,8),temp(:,12),];
        
        e_xx=G(1,:);
        e_yy=G(5,:);
        e_zz=G(9,:);
        
        e_xy=G(2,:)+G(4,:);
        e_xz=G(3,:)+G(7,:);
        e_yz=G(6,:)+G(8,:);
        
        shear_vector=[e_xx;e_yy;e_zz;e_xy;e_xz;e_yz];
        
        A=zeros(12);
        syms('E','v','x_','y_','real');           
        
        C_inv=blkdiag(-(my_value)*ones(3)+(my_value+1)*eye(3),my_frac*eye(3))/E_value;
        C=inv(C_inv)*det(C_inv)        
     
        for i=1:12
            for j=1:12
                A(i,j)=abs(det_J_t)*dot(shear_vector(:,i),C*shear_vector(:,j));
            end
        end
        
        
        Mass_mat=zeros(12);
                
        for i=1:4
                RHS_pre_phi_2=@(x,y,z) ([x,y,z]-x_4)*G(1:3,1:3)*basis_f(:,i)+corr_vec(i);
            for j=1:4
                RHS_pre_phi_1 =@(x,y,z) ([x,y,z]-x_4)*G(1:3,1:3)*basis_f(:,j)+corr_vec(j);
                Mass_mat((i-1)*3+1,(j-1)*3+1)=gauss_quad_3(x_1,x_2,x_3,x_4,5,@(x,y,z) RHS_pre_phi_1(x,y,z)*RHS_pre_phi_2(x,y,z));
                Mass_mat((i-1)*3+2,(j-1)*3+2)=Mass_mat((i-1)*3+1,(j-1)*3+1);
                Mass_mat((i-1)*3+3,(j-1)*3+3)=Mass_mat((i-1)*3+1,(j-1)*3+1);
            end
        end
        Mass_mat=density*Mass_mat
        
        
        
        
end


