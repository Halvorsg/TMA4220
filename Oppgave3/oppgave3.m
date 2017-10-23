function [u_sol,u_exact] = oppgave3(N)

addpath ..\Grids
addpath ..\Oppgave1
%% Initial conditions
theta = 3*pi/2;
%% Get triangle
[p,tri,edges] = getSlice(N,theta);
TR = triangulation(tri, p);
% triplot(TR)
% hold on
%% Functions
f = @(x,y) 16*pi^2*x.*y.*(x.^2+y.^2).*sin(2*pi*(x.^2+y.^2))-24*x.*y.*pi.*cos(2*pi*(x.^2+y.^2));
exact_solution = @(x,y) x.*y.*sin(2*pi.*(x.^2+y.^2));
gx = @(x,y) -x.*sin(2*pi*x.^2);
gy = @(x,y)  y.*sin(2*pi*y.^2);
%% Pre-allocating
A = spalloc(length(p),length(p),10*length(p));
F = zeros(length(p),1);
phi = [eye(2),[-1;-1]];
%% Taking care of Neumann boundary
p(abs(p)<100*eps) = 0;      % Setting all values of p approx 0 to 0
edges1 = p(edges(:,1),:);   % Position of first node in edge pair
edges2 = p(edges(:,2),:);   % Position of second node in edge pair
Neumann_edges = zeros(2,2); % Pre-allocating table of Neumann edges
% i = 1;
for i = 1:length(edges)
% plot([edges1(i,1),edges2(i,1)],[edges1(i,2),edges2(i,2)],'r')
    switch is_Neumann_edge(edges1(i,:),edges2(i,:))
        case 1
            %Vertical <=> x = 0
            Fn = quadratureLine2D(edges1(i,:),edges2(i,:),3,gy);            % Line integral for Neumann boundary
            F(edges(i,:)) = F(edges(i,:)) + Fn;                             % Adding to load vector
%             plot([edges1(i,1),edges2(i,1)],[edges1(i,2),edges2(i,2)],'g')   % Plotting edge in triplot for visualization
            Neumann_edges(i,:) = edges(i,:);                                % Adding edge(i) to list of Neumann edges
            
        case 2
            %Horizontal <=> y = 0
            Fn = quadratureLine2D(edges1(i,:),edges2(i,:),3,gx);            % Line integral for Neumann boundary
            F(edges(i,:)) = F(edges(i,:)) + Fn;                             % Adding to load vector
%             plot([edges1(i,1),edges2(i,1)],[edges1(i,2),edges2(i,2)],'g')   % Plotting edge in triplot for visualization
            Neumann_edges(i,:) = edges(i,:);                                % Adding edge(i) to list of Neumann edges
        case 0
            %do nothing
    end
end
% x = ismember(edges,Neumann_edges);
% edges = edges(sum(x,2) < 2 , : );
edges = edges(sum(ismember(edges,Neumann_edges),2) < 2 , : );               % Removing Neumann edges from Dirchlet edges
%% Creating A and finishing F
for i = 1:length(tri)
    p1 = p(tri(i,1),:)'; p2 = p(tri(i,2),:)'; p3 = p(tri(i,3),:)'; 
    K = [p1-p3 , p2-p3];
    jacDet = abs(det(K));
%Stiffness element matrix - Ak
    G = K'\phi;
    Ak = jacDet*(G'*G)/2;
%Load elemen vector - Fk
    Fk = quadrature2D(p1,p2,p3,3,f);      
%Adding into Stiffnessmatrix and load vector
    A(tri(i,:),tri(i,:)) = A(tri(i,:),tri(i,:)) + Ak;
    F(tri(i,:)) = F(tri(i,:)) + Fk;
end
%% Adding and removing boundary
u_sol = zeros(size(p(:,1)));

y = ismember(tri,edges);
z = sum(y,2);
x = tri(z<=1,:);
m = ~ismember(x,edges);
n = x(m);
inner_vertices = unique(n);
A = A(inner_vertices,inner_vertices);
F = F(inner_vertices);

%% Solving problem and plotting
u = A\F;
u_sol(inner_vertices) = u;
u_exact = exact_solution(p(:,1),p(:,2));
figure
trimesh(tri,p(:,1),p(:,2),u_sol)
s = sprintf('Approximate solution with N = %i', N);
title(s)
figure
trimesh(tri,p(:,1),p(:,2),u_exact)
title('Exact solution')
end
    
    
    
    
    