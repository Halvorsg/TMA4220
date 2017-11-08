function [ int_value ] = gauss_quad_2(point_1, point_2, point_3, nr_of_int_points, func_to_int )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


        switch nr_of_int_points
            case 1
              lambda = [1/3; 1/3; 1/3];
              weight = 1;
            case 3
              lambda = [1/2 1/2 0;1 0 2; 0 2 2];
              weight = [1/3 1/3 1/3];
            case 4
              lambda = [1/3 3/5 1/5 1/5;1/3 1/5 3/5 1/5;1/3 1/5 1/5 3/5];
              weight = [-9/16 25/48 25/48 25/48];
            otherwise
                fprintf("invalod choice of int points");
                error('#1')
        end   
        
       
       
        points_eval=[transpose(point_1) transpose(point_2) transpose(point_3)]*lambda;
        int_value=dot(arrayfun(func_to_int,points_eval(1,:),points_eval(2,:)),weight);
        
        J_t=[transpose(point_1-point_3), transpose(point_2-point_3)];
        det_J_t=det(J_t);
  
        int_value=int_value*det_J_t;

end

