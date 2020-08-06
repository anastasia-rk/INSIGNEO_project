clc; clear; local_init;
%% Testing analytical expressions
x = [0:0.1:4];
y0 = bsplinederiv(4,0,4,x,1);
figure; plot(x,y0);
% y1 = bsplinederiv(4,0,4,x,2);
% figure; plot(x,y1);
% y3 = bsplinederiv(4,0,4,x,3);
% figure; plot(x,y3);
% y4 = bsplinederiv(4,0,4,x,4);
% figure; plot(x,y4);

% y = [0:0.2:4];
% dfdx = partial_deriv(4,0,4,0,4,x,y,1);
% figure; plot3(x,y,dfdx);
% partial derivative w.r.t to x
for i=1:length(x)
    y(i,:) = [0:0.1:4];
    x_plot(:,i) = x(i);
    dfdx(i,:) = partial_deriv(4,0,4,0,4,x(i),y(i,:),1);
end
figure; surf(x_plot',y',dfdx');
colormap(my_map);
xlabel('x'); ylabel('y');zlabel('df/dx')
clear y x
% partial derivative w.r.t. to y
y = [0:0.1:4];
for i=1:length(y)
    x(i,:) = [0:0.1:4];
    y_plot(:,i) = y(i);
    dfdy(i,:) = partial_deriv(4,0,4,0,4,y(i),x(i,:),1);
end
figure; surf(x,y_plot,dfdy);
colormap(my_map);
xlabel('x'); ylabel('y');zlabel('df/dy')
%%
% patial derivatives
function [dfdx] = partial_deriv(order,ax,bx,ay,by,xval,yval,dorder)
% x - coordinate w.r.t. which we take partial derivative
dbx = bsplinederiv(order,ax,bx,xval,dorder);
by  = bsplinexval(order,ay,by,yval);
dfdx = dbx.*by;
end

function bx = bsplinexval (order,a,b,xval)
x = linspace(a,b,order+1);
%get bspline in ppform
pp = bspline(x);
%bx is the value of the spline at coordinate xval 
index = find((xval > a) & (xval < b));
bx = zeros(size(xval));
bx(index) = fnval(pp,xval(index));
end 

function bx = bsplinederiv(order,a,b,xval,dorder)
x = linspace(a,b,order+1);
%get bspline in ppform
pp = bspline(x);
%derivative of order with order of derivative 
db1 = fnder(pp,dorder);
% evaluate the derivative at xval
bx = fnval(db1,xval);
end 