function[grad] = gradient_bspline(position,knots,order)
% Gradient of tensor product b-spline grid evaluated at the single point
% position (2x1) - corridnates of the point of interest
% knots (2Lx1)   - support borders for b-splines
% order          - b-spline order
xval = position(1,:);
yval = position(2,:);
iGrad = 0;
for i=1:2:length(knots)-1
    ax = knots(1,i);
    bx = knots(1,i+1);
    ay = knots(2,i);
    by = knots(2,i+1);
    dfdx = partial_der(order,ax,bx,ay,by,xval,yval);
    dfdy = partial_der(order,ay,by,ax,bx,yval,xval);
    iGrad = iGrad + 1;
    grad(:,iGrad) = [dfdx;dfdy];
end