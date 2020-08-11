function bx = bsplinexval (order,a,b,xval)
x = linspace(a,b,order+1);
%get bspline in ppform
pp = bspline(x);
%bx is the value of the spline at coordinate xval 
index = find((xval > a) & (xval < b));
bx = zeros(size(xval));
bx(index) = fnval(pp,xval(index));
end 