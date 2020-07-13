%% Testing for the modified 2d b-spline 
clc; clear;
by = 10; bx = 10;
xval1 = [0:0.2:20];
yval1 = [0:0.2:20];
[xx,yy] = meshgrid(xval1,yval1);
zz1 = tensorproductbspline(4,0,bx,0,by,xx,yy);
xval2 = [0:0.2:20];
yval2 = [0:0.2:20];
zz2 = tensorproductbspline(4,4,bx+4,4,by+4,xx,yy);
%Plot surface 
figure;
surf(xx,yy,zz1)
hold on;
surf(xx,yy,zz2)
alpha(0.5)
%Labels
xlabel('x');
ylabel('y');
zlabel('z');
%% Functions
function zz = tensorproductbspline(order,ax,bx,ay,by,xval,yval)
% Calling previous function to form bsplines from start a to end b, with
% specified order and value xval
% bxx and byy plot the surface
% bxx = bspline123(order,ax,bx);
% byy = bspline123(order,ay,by);
%pp and tt plot the point of interest 
pp = bsplinexval(order,ax,bx,xval);
tt = bsplinexval(order,ay,by,yval);
%Creating meshgrid of bxx and byy
% [xx,yy] = meshgrid(bxx,byy);
%Define relationship between x, y and z
zz = (pp).*(tt);
end 
%knots from start point a to end point b o is order of the spline 
%and xval is the x coordinate to evluate at
function bx = bsplinexval (order,a,b,xval)
x = linspace(a,b,order+1);
%get bspline in ppform
pp = bspline(x);
%bx is the value of the spline at coordinate xval 
index = find((xval > a) & (xval < b))
bx = zeros(size(xval));
bx(index) = fnval(pp,xval(index));
end 
% out put is bx value for xval