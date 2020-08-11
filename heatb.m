function heatb(order,ax,bx,ay,by,xval,yval,theta)
%make meshgrid of x and y values 
[xx,yy] = meshgrid(xval,yval);
%call tensor product bspline function to get zz- the tenosr prodcut of x
%and y 
zz = tensorproductbspline(order,ax,bx,ay,by,xx,yy);
Z_plot = zz*theta;
% plotting heat map
Z_min = min(Z_plot);
Z_min_min = min(Z_min);
a_min = min(min(Z_plot));
a_max = max(max(Z_plot));
p = pcolor(xx,yy,Z_plot);
view(2);
shading interp
caxis manual
alpha(p,0.8);
caxis([a_min a_max]);
% plot the surface 
surf(xx,yy,Z_plot)
%Labels
title('Heat plot');
colorbar; %('north','Color','k','FontSize',12);
xlabel('x');
ylabel('y');
zlabel('z');
end 
