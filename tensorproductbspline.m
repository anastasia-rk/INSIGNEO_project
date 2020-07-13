
function zz = tensorproductbspline(order,ax,bx,ay,by,xxval,yyval)
% Calling previous function to form bsplines from start a to end b, with
% specified order and value xval
%pp and tt plot the point of interest 
pp = bsplinexval(order,ax,bx,xxval);
tt = bsplinexval(order,ay,by,yyval);
%Creating meshgrid of bxx and byy
%Define relationship between x, y and z
zz = (pp).*(tt);
end 
