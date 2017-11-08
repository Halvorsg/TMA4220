function [ u_sol,u_exact] = temperature(N)
addpath('..\Grids');
addpath('..\Oppgave1');
%% Get triangle
[p,tri,edges] = getPlate(N);
%% Functions
u1 = @(x,y) x;
u2 = @(x,y) 10;
f = @(x,y) 2*(x^2 - 1)*(y^2 + 10*y - 2);                              % Right hand side
gtop = @(x,y) 2*(x^2-1);                            % Neumann top
gbottom = @(x,y)  2*(x^2-1);                        % Neumann bottom % Be carefull with the normal vector
exact_solution = @(x,y) (x.^2 - 1).*(y.^2 - 1); % Exact solution
%% Pre-allocating
A = spalloc(length(p),length(p),10*length(p));
M = spalloc(length(p),length(p),10*length(p));
F = zeros(length(p),1);
phi = [eye(2),[-1;-1]];
% Taking care of Neumann boundary
p(abs(p)<100*eps) = 0;      % Setting all values of p approx 0 to 0
edges1 = p(edges(:,1),:);   % Position of first node in edge pair
edges2 = p(edges(:,2),:);   % Position of second node in edge pair
Neumann_edges = zeros(2*N,2); % Pre-allocating table of Neumann edges
% %% Neumane edges
% 
% % i = 1;
% for i = 1:length(edges)
% % plot([edges1(i,1),edges2(i,1)],[edges1(i,2),edges2(i,2)],'r')
%     p1 = edges1(i,:); p2 = edges2(i,:);
%     if p1(2) == -1 && p2(2) == -1
%             %Vertical <=> x = 0
%             Fn = quadratureLine2D(edges1(i,:),edges2(i,:),3,gbottom);            % Line integral for Neumann boundary
%             F(edges(i,:)) = F(edges(i,:)) + Fn;                             % Adding to load vector
% %             plot([edges1(i,1),edges2(i,1)],[edges1(i,2),edges2(i,2)],'g')   % Plotting edge in triplot for visualization
%             Neumann_edges(i,:) = edges(i,:);                                % Adding edge(i) to list of Neumann edges
%             
%     elseif p1(2) == 1 && p2(2) == 1
%             %Horizontal <=> y = 0
%             Fn = quadratureLine2D(edges1(i,:),edges2(i,:),3,gtop);            % Line integral for Neumann boundary
%             F(edges(i,:)) = F(edges(i,:)) + Fn;                             % Adding to load vector
% %             plot([edges1(i,1),edges2(i,1)],[edges1(i,2),edges2(i,2)],'g')   % Plotting edge in triplot for visualization
%             Neumann_edges(i,:) = edges(i,:);                                % Adding edge(i) to list of Neumann edges
%             %do nothing
%     end
% end
% 
% edges = edges(sum(ismember(edges,Neumann_edges),2) < 2 , : );               % Removing Neumann edges from Dirchlet edges
%% Creating A and finishing F
for i = 1:length(tri)
    p1 = p(tri(i,1),:)'; p2 = p(tri(i,2),:)'; p3 = p(tri(i,3),:)'; 
    K = [p1-p3 , p2-p3];
    jacDet = abs(det(K));
%Stiffness element matrix - Ak
    G = K'\phi;
    Ak = jacDet*(G'*G)/2;
% Mass matrix
 
    
    Mk=zeros(3);
    cor_vec=[0 0 1];
    for it=1:3;
        for j=1:3;
            func_to_int = @(x,y)dot([u1(x,y),u2(x,y)],(([x,y]-p3')*G(1:2,1:2))'.*G(:,it))*(([x,y]-p3')*G(1:2,1:2)*phi(:,j)+cor_vec(j));
            Mk(it,j)=gauss_quad_2(p1',p2',p3',4,func_to_int);
        end
    end
    
   % phi_u1 = quadrature2D(p_1,p_2,p_3,3,u1);
   % phi_u2 = quadrature2D(p_1,p_2,p_3,3,u2);
   % phi_u = [phi_u1,phi_u2];
   % Mk = jacDet*phi_u*G;
%Load elemen vector - Fk
    Fk = quadrature2D(p1,p2,p3,3,f);      
%Adding into Stiffnessmatrix and load vector
    A(tri(i,:),tri(i,:)) = A(tri(i,:),tri(i,:)) + Ak;
    M(tri(i,:),tri(i,:)) = M(tri(i,:),tri(i,:)) + Mk;
    F(tri(i,:)) = F(tri(i,:)) + Fk;
end
%% Adding and removing boundary
u_sol = zeros(size(p(:,1)));

y = (~ismember(tri,edges)).*tri;
inner_vertices = unique(y);
inner_vertices = inner_vertices(2:end);

%% Lifting and solving
edges_hot = edges(p(edges) == -1);
edges_cold = edges(p(edges) == 1);

gr = zeros(length(p),1);
gr(edges_hot) = 0; % Kelvin
gr(edges_cold) = 0; % Kelvin

G = F- (A+M)*gr;
A = A(inner_vertices,inner_vertices);
M = M(inner_vertices,inner_vertices);
F = G(inner_vertices);
u = (A+M)\F;

u_sol(inner_vertices) = u;
u_sol = u_sol+gr;
u_exact = exact_solution(p(:,1),p(:,2));

%% Plotting
figure
trimesh(tri,p(:,1),p(:,2),u_sol)
s = sprintf('Approximate solution with N = %i', N);
xlabel('x')
ylabel('y')
title(s)

figure
trimesh(tri,p(:,1),p(:,2),u_exact)
title('Exact solution')

figure
trimesh(tri,p(:,1),p(:,2),u_exact-u_sol)
xlabel('x')
ylabel('y')
end

