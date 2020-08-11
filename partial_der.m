function [dfdx] = partial_der(order,ax,bx,ay,by,xval,yval)
% x - coordinate w.r.t. which we take partial derivative
dbx = bsplinederiv(order,ax,bx,xval,1);
by  = bsplinexval(order,ay,by,yval);
% expression for first order derivative only
dfdx = dbx.*by;
end