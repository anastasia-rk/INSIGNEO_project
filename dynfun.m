function[out] = dynfun(x,A,B,theta,knots,order)
beta = gradient_bspline(x(1:2),knots,order);
out = A*x + B*beta*theta;