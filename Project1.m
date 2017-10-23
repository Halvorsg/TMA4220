function Project1()
% Run this function to se how each task is excecuted
% I would recommend to run each section by itself
addpath Oppgave1
addpath Oppgave2
addpath Oppgave3
addpath Oppgave4
addpath Grids
%% Task 1

g = @(x) exp(x);
I1 = quadrature1D(1,2,4,g);
g = @(x,y) log(x+y);
p1 = [1;0]; p2 = [3;1]; p3 = [3;2];
I2 = quadrature2D_Show(p1,p2,p3,4,g);
g = @(x,y,z) exp(x);
p1 = [0,0,0]'; p2 = [0,2,0]'; p3 = [0,0,2]'; p4 = [2,0,0]';
I3 = quadrature3D_Show(p1,p2,p3,p4,4,g);
%% Task 2
N = 1000; % Number of points
[approx,exact] = oppgave2(N);
err2 = 1/sqrt(length(approx))*norm(approx-exact);
errInf = norm(approx-exact,'inf');
fprintf('Error in Task 2 for N = %i \n L2-Norm of error = %f\n infinity-norm of error = %f\n\n',length(approx),err2,errInf)
%% Task 3
N = 1000; %Number of points
[approx,exact] = oppgave3(N);
err2 = 1/sqrt(length(approx))*norm(approx-exact);
errInf = norm(approx-exact,'inf');
fprintf('Error in Task 3 for N = %i \n L2-Norm of error = %f\n infinity-norm of error = %f\n\n',length(approx),err2,errInf)
%% Task 4a
[approx,exact] = oppgave4(10);
err2 = 1/sqrt(length(approx))*norm(approx-exact);
errInf = norm(approx-exact,'inf');
fprintf('Error in Task 4a for N = %i \n L2-Norm of error = %f\n infinity-norm of error = %f\n\n',length(approx),err2,errInf)
%% Task 4b and 4c
[approx,exact] = oppgave4b(1000);
err2 = 1/sqrt(length(approx))*norm(approx-exact);
errInf = norm(approx-exact,'inf');
fprintf('Error in Task 4b for N = %i \n L2-Norm of error = %f\n infinity-norm of error = %f\n\n',length(approx),err2,errInf)
end