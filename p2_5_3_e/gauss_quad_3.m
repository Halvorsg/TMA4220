function [ int_value ] = gauss_quad_3(point_1, point_2, point_3, point_4, nr_of_int_points, func_to_int)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    switch nr_of_int_points
                case 1
                  lambda = [1/4; 1/4; 1/4; 1/4];
                  weight = 1;
                case 4
                  lambda = [0.5854102 0.1381966 0.1381966 0.1381966; 0.1381966 0.5854102  0.1381966  0.1381966; 0.1381966  0.1381966 0.5854102  0.1381966;  0.1381966  0.1381966  0.1381966 0.5854102];
                  weight = [1/4 1/4 1/4 1/4];
                case 5
                  lambda = [1/4 1/2 1/6 1/6 1/6; 1/4 1/6 1/2 1/6 1/6; 1/4 1/6 1/6 1/2 1/6; 1/4 1/6 1/6 1/6 1/2];
                  weight = [-4/5 9/20 9/20 9/20 9/20];
                otherwise
                    fprintf("invalod choice of int points");
                    error('#1')
    end

points_eval=[transpose(point_1) transpose(point_2) transpose(point_3) transpose(point_4)]*lambda;
        int_value=dot(arrayfun(func_to_int,points_eval(1,:), points_eval(2,:), points_eval(3,:)),weight);

        J_t=[transpose(point_1-point_4), transpose(point_2-point_4), transpose(point_3-point_4)];
        det_J_t=det(J_t);
  
        int_value=int_value*(1/3)*(abs(det_J_t));
        
end

