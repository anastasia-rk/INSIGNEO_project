function tensorproductbspline(order,ax,bx,ay,by,xval,yval)
% Calling previous function to form bsplines from start a to end b, with
% specified order and value xval
% bxx and byy plot the surface
bxx = bspline123(order,ax,bx);
byy = bspline123(order,ay,by);
%pp and tt plot the point of interest 
pp = bsplinexval(order,ax,bx,xval);
tt = bsplinexval(order,ay,by,yval);
%Creating meshgrid of bxx and byy
[xx,yy] = meshgrid(bxx,byy);
%Define relationship between x, y and z
z = (xx).*(yy);
zz = (pp).*(tt)
%Plot surface 
figure;
surf(z)
hold on 
plot3(xval,yval,zz,'*r')
hold off
%Labels
xlabel('x');
ylabel('y');
zlabel('z');
end 
%knots from start point a to end point b o is order of the spline 
%and xval is the x coordinate to evluate at
function bx = bspline123 (order,a,b)
x = linspace(a,b,order+1);
%get bspline in ppform
pp = bspline(x);
%bx is the value of the spline at coordinate xval 
bx = fnval(pp,x)
end 
% out put is bx value for xval

%knots from start point a to end point b o is order of the spline 
%and xval is the x coordinate to evluate at
function bx = bsplinexval (order,a,b,xval)
x = linspace(a,b,order+1);
%get bspline in ppform
pp = bspline(x);
%bx is the value of the spline at coordinate xval 
bx = fnval(pp,xval)
end 
% out put is bx value for xval
