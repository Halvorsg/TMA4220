function contour_plot_100100( V,real_inner_vertices,p ,k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
scaling = 10;

x=zeros(length(p),1);
y=zeros(length(p),1);
z=zeros(length(p),1);

step=1;
  for i=1:length(real_inner_vertices)
    if p(real_inner_vertices(i),3)>0.48
    x(step)=p(real_inner_vertices(i),1)+V(i*3-2,k)*scaling;
    y(step)=p(real_inner_vertices(i),2)+V(i*3-1,k)*scaling;
    z(step)=p(real_inner_vertices(i),3)+V(i*3,k)*scaling;
    step=step+1;
    end
  end
    x=x(1:step-1);
    y=y(1:step-1);
    z=z(1:step-1);
   
[m_x,m_y] = meshgrid(0:1:100, 0:1:100);
test_plot = griddata(x,y,z,m_x,m_y);%   
contour(m_x,m_y,test_plot);
Animation=figure(1);

end

