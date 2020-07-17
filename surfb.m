function surfb(order,ax,bx,ay,by,xval,yval,theta)
%make meshgrid of x and y values 
[xx,yy] = meshgrid(xval,yval);
%call tensor product bspline function to get zz- the tenosr prodcut of x
%and y 
zz = tensorproductbspline(order,ax,bx,ay,by,xx,yy);
ZZ = zz*theta;
% plot the surface 
surf(xx,yy,ZZ)
%Labels
xlabel('x');
ylabel('y');
zlabel('z');
end 
