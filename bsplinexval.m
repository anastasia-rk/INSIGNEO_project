%knots from start point a to end point b o is order of the spline 
%and xval is the x coordinate to evluate at
function bx = bsplinexval (order,a,b,xval)
x = linspace(a,b,order+1);
%get bspline in ppform
pp = bspline(x);
%bx is the value of the spline at coordinate xval 
%find xval bigger than a and smaller than b
index = find((xval > a) & (xval < b))
%creates zero matrix the size of xval
bx = zeros(size(xval));
%inserts value between a and b leaving the rest at zero
bx(index) = fnval(pp,xval(index));
end 