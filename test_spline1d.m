%% Testing for the modified 1d b-spline 
x = [0:0.2:7];
y = bsplinexval(4,0,4,x);
figure; plot(x,y);

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