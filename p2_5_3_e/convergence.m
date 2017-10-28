function convergence()
N = [5,10,15];
e = zeros(length(N),3);

for i = 1:length(N)
    tic
    [ux,uy,uz,u_exact] = task53e(N(i));
    e(i,:) = [max(abs(ux-u_exact)),max(abs(uy-u_exact)),max(abs(uz-u_exact))];
    toc
end

loglog(N,e,'-*')
save('C:\Users\halvo\Documents\MATLAB\TMA4220\TMA4220\p2_5_3_e\errors.mat','e')
save('C:\Users\halvo\Documents\MATLAB\TMA4220\TMA4220\p2_5_3_e\gridSize.mat','N')

end