%knots from start point a to end point b o is order of the spline 
%and xval is the x coordinate to evluate at
function bx = bsplinexval(order,a,b,xval,dorder)
x = linspace(a,b,order+1);
%get bspline in ppform
pp = bspline(x);
%derivative of order with order of derivative 
db1 = fnder(pp,dorder)
p11 = fnplt(db1,'j')
%plot p11 separate rows
plot(p11(1,:),p11(2,:))
end 