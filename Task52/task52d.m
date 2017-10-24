function [u_sol,u_exact] = task52d(N)

addpath ..\Grids
addpath ..\Oppgave1
%% Get triangulation
[p,tri,edge] = getPlate(N);

%% Set initial values and make functions for Ak and Fk
E = 1;
v = 0;
fx = @(x,y) E/(1-v^2) * [-2*y.^2 + x.^2*(v-1) - 2*x.*y*(v+1) + 3 - v];
fy = @(x,y) E/(1-v^2) * [-2*x.^2 + y.^2*(v-1) - 2*x.*y*(v+1) + 3 - v];
Ak =  @(E,v,p1,p2,p3)reshape([(E.*(p2(1).*p3(1).*-2.0-p2(2).*p3(2).*4.0-v.*p2(1).^2-v.*p3(1).^2+p2(1).^2+p3(1).^2+p2(2).^2.*2.0+p3(2).^2.*2.0+v.*p2(1).*p3(1).*2.0).*(-1.0./4.0))./((v.^2-1.0).*(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2))),(E.*(p2(1)-p3(1)).*(p2(2)-p3(2)).*(1.0./4.0))./((v-1.0).*(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2))),(E.*((p1(2)-p3(2)).*(p2(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-(v.*(1.0./2.0)-1.0./2.0).*(p1(1)-p3(1)).*(p2(1)-p3(1)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*((v.*(1.0./2.0)-1.0./2.0).*(p2(1)-p3(1)).*(p1(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-v.*(p1(1)-p3(1)).*(p2(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),-(E.*((p1(2)-p2(2)).*(p2(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-(v.*(1.0./2.0)-1.0./2.0).*(p1(1)-p2(1)).*(p2(1)-p3(1)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),-(E.*((v.*(1.0./2.0)-1.0./2.0).*(p2(1)-p3(1)).*(p1(2)-p2(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-v.*(p1(1)-p2(1)).*(p2(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*(p2(1)-p3(1)).*(p2(2)-p3(2)).*(1.0./4.0))./((v-1.0).*(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2))),(E.*(p2(1).*p3(1).*-4.0-p2(2).*p3(2).*2.0-v.*p2(2).^2-v.*p3(2).^2+p2(1).^2.*2.0+p3(1).^2.*2.0+p2(2).^2+p3(2).^2+v.*p2(2).*p3(2).*2.0).*(-1.0./4.0))./((v.^2-1.0).*(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2))),(E.*((v.*(1.0./2.0)-1.0./2.0).*(p1(1)-p3(1)).*(p2(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-v.*(p2(1)-p3(1)).*(p1(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*((p1(1)-p3(1)).*(p2(1)-p3(1)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-(v.*(1.0./2.0)-1.0./2.0).*(p1(2)-p3(2)).*(p2(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),-(E.*((v.*(1.0./2.0)-1.0./2.0).*(p1(1)-p2(1)).*(p2(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-v.*(p2(1)-p3(1)).*(p1(2)-p2(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),-(E.*((p1(1)-p2(1)).*(p2(1)-p3(1)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-(v.*(1.0./2.0)-1.0./2.0).*(p1(2)-p2(2)).*(p2(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*((p1(2)-p3(2)).*(p2(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-(v.*(1.0./2.0)-1.0./2.0).*(p1(1)-p3(1)).*(p2(1)-p3(1)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*((v.*(1.0./2.0)-1.0./2.0).*(p1(1)-p3(1)).*(p2(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-v.*(p2(1)-p3(1)).*(p1(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*(p1(1).*p3(1).*-2.0-p1(2).*p3(2).*4.0-v.*p1(1).^2-v.*p3(1).^2+p1(1).^2+p3(1).^2+p1(2).^2.*2.0+p3(2).^2.*2.0+v.*p1(1).*p3(1).*2.0).*(-1.0./4.0))./((v.^2-1.0).*(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2))),(E.*(p1(1)-p3(1)).*(p1(2)-p3(2)).*(1.0./4.0))./((v-1.0).*(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2))),(E.*((p1(2)-p2(2)).*(p1(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-(v.*(1.0./2.0)-1.0./2.0).*(p1(1)-p2(1)).*(p1(1)-p3(1)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*((v.*(1.0./2.0)-1.0./2.0).*(p1(1)-p3(1)).*(p1(2)-p2(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-v.*(p1(1)-p2(1)).*(p1(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*((v.*(1.0./2.0)-1.0./2.0).*(p2(1)-p3(1)).*(p1(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-v.*(p1(1)-p3(1)).*(p2(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*((p1(1)-p3(1)).*(p2(1)-p3(1)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-(v.*(1.0./2.0)-1.0./2.0).*(p1(2)-p3(2)).*(p2(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*(p1(1)-p3(1)).*(p1(2)-p3(2)).*(1.0./4.0))./((v-1.0).*(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2))),(E.*(p1(1).*p3(1).*-4.0-p1(2).*p3(2).*2.0-v.*p1(2).^2-v.*p3(2).^2+p1(1).^2.*2.0+p3(1).^2.*2.0+p1(2).^2+p3(2).^2+v.*p1(2).*p3(2).*2.0).*(-1.0./4.0))./((v.^2-1.0).*(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2))),(E.*((v.*(1.0./2.0)-1.0./2.0).*(p1(1)-p2(1)).*(p1(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-v.*(p1(1)-p3(1)).*(p1(2)-p2(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*((p1(1)-p2(1)).*(p1(1)-p3(1)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-(v.*(1.0./2.0)-1.0./2.0).*(p1(2)-p2(2)).*(p1(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),-(E.*((p1(2)-p2(2)).*(p2(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-(v.*(1.0./2.0)-1.0./2.0).*(p1(1)-p2(1)).*(p2(1)-p3(1)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),-(E.*((v.*(1.0./2.0)-1.0./2.0).*(p1(1)-p2(1)).*(p2(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-v.*(p2(1)-p3(1)).*(p1(2)-p2(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*((p1(2)-p2(2)).*(p1(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-(v.*(1.0./2.0)-1.0./2.0).*(p1(1)-p2(1)).*(p1(1)-p3(1)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*((v.*(1.0./2.0)-1.0./2.0).*(p1(1)-p2(1)).*(p1(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-v.*(p1(1)-p3(1)).*(p1(2)-p2(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*(p1(1).*p2(1).*-2.0-p1(2).*p2(2).*4.0-v.*p1(1).^2-v.*p2(1).^2+p1(1).^2+p2(1).^2+p1(2).^2.*2.0+p2(2).^2.*2.0+v.*p1(1).*p2(1).*2.0).*(-1.0./4.0))./((v.^2-1.0).*(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2))),(E.*(p1(1)-p2(1)).*(p1(2)-p2(2)).*(1.0./4.0))./((v-1.0).*(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2))),-(E.*((v.*(1.0./2.0)-1.0./2.0).*(p2(1)-p3(1)).*(p1(2)-p2(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-v.*(p1(1)-p2(1)).*(p2(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),-(E.*((p1(1)-p2(1)).*(p2(1)-p3(1)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-(v.*(1.0./2.0)-1.0./2.0).*(p1(2)-p2(2)).*(p2(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*((v.*(1.0./2.0)-1.0./2.0).*(p1(1)-p3(1)).*(p1(2)-p2(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-v.*(p1(1)-p2(1)).*(p1(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*((p1(1)-p2(1)).*(p1(1)-p3(1)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))-(v.*(1.0./2.0)-1.0./2.0).*(p1(2)-p2(2)).*(p1(2)-p3(2)).*(p1(1).*p2(2).*(1.0./2.0)-p2(1).*p1(2).*(1.0./2.0)-p1(1).*p3(2).*(1.0./2.0)+p3(1).*p1(2).*(1.0./2.0)+p2(1).*p3(2).*(1.0./2.0)-p3(1).*p2(2).*(1.0./2.0))).*1.0./(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)).^2)./(v.^2-1.0),(E.*(p1(1)-p2(1)).*(p1(2)-p2(2)).*(1.0./4.0))./((v-1.0).*(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2))),(E.*(p1(1).*p2(1).*-4.0-p1(2).*p2(2).*2.0-v.*p1(2).^2-v.*p2(2).^2+p1(1).^2.*2.0+p2(1).^2.*2.0+p1(2).^2+p2(2).^2+v.*p1(2).*p2(2).*2.0).*(-1.0./4.0))./((v.^2-1.0).*(p1(1).*p2(2)-p2(1).*p1(2)-p1(1).*p3(2)+p3(1).*p1(2)+p2(1).*p3(2)-p3(1).*p2(2)))],[6,6]);

%% Preallocate and iterate over all elements
N = length(tri);
A = spalloc(2*N,2*N,10*N);
F = zeros(2*N,1);


for i = 1:length(tri)
    p1 = p(tri(i,1),:)'; p2 = p(tri(i,2),:)'; p3 = p(tri(i,3),:)';
    %Finding the right indices for A by i = 2*(�-1) + d 
    tempIndices = [2*(tri(i,:)-1)+1;2*(tri(i,:)-1)+2];
    indices = tempIndices(:);
    A(indices,indices) = A(indices,indices) + Ak(E,v,p1,p2,p3);
    % Finding the quadrature for fx and fy separate, easier to use existing
    % code this way
    Fkx = quadrature2D(p1,p2,p3,3,fx);
    Fky = quadrature2D(p1,p2,p3,3,fy);
    Fk = [Fkx';Fky'];
    Fk = Fk(:);
    F(indices) = F(indices) + Fk;
end
%% Implementing boundary conditions and solving
y = (~ismember(tri,edge)).*tri;
inner_vertices_2 = unique(y);
inner_vertices_2 = inner_vertices_2(2:end);

tempIV = [2*(inner_vertices_2-1)+1;2*(inner_vertices_2-1)+2];
inner_vertices = tempIV(:);

A = A(inner_vertices,inner_vertices);
F = F(inner_vertices);


u = A\F;
u_sol = zeros(2*length(p),1);
u_sol(inner_vertices) = u;
ux = u_sol(1:2:end);
uy = u_sol(2:2:end);

%% Plotting
figure
trimesh(tri,p(:,1),p(:,2),ux)
figure
trimesh(tri,p(:,1),p(:,2),uy)
exact_solutionx = @(x,y) (x.^2-1).*(y.^2-1);
ux_exact = exact_solutionx(p(:,1),p(:,2));
figure
trimesh(tri,p(:,1),p(:,2),ux_exact)
end
    