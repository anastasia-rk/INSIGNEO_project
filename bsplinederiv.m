function bx = bsplinederiv(order,a,b,xval,dorder)
x = linspace(a,b,order+1);
%get bspline in ppform
pp = bspline(x);
%derivative of order with order of derivative 
db1 = fnder(pp,dorder);
% evaluate the derivative at xval
bx = fnval(db1,xval);
% p11 = fnplt(db1,'j');
% 
% %plot p11 separate rows
% figure;
% plot(p11(1,:),p11(2,:))
end 