
% patial derivatives
function [dfdx] = partial_derivative(order,ax,bx,ay,by,xval,yval,dorder)
% x - coordinate w.r.t. which we take partial derivative
dbx = bsplinederiv(order,ax,bx,xval,dorder);
by  = bsplinexval(order,ay,by,yval);
dfdx = dbx.*by;
end