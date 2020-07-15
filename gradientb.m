function gradientb(order,ax,bx,ay,by,xval,yval,theta)
%make meshgrid of x and y values 
[xx,yy] = meshgrid(xval,yval);
%call tensor product bspline function to get zz- the tenosr prodcut of x
%and y 
zz = tensorproductbspline(order,ax,bx,ay,by,xx,yy);
Z_plot = zz*theta;
[dx, dy] = gradient(Z_plot, 10, 10);
Z_min = min(Z_plot);
Z_min_min = min(Z_min);
contour(xx, yy, Z_plot), hold on    % plot level surfaces
quiver(xx, yy, dx, dy, 1, 'k'); % plot vector field
% plot the surface 
title('Gradient plot');
%Labels
xlabel('x');
ylabel('y');
zlabel('z');
end 
