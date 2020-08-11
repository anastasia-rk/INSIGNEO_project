x = [0:0.1:4];
y0 = bsplinederiv(4,0,4,x,0)
plot(x,y0,'.-');
hold on 
y  = bsplinederiv(4,0,4,x,1)
plot(x,y);
hold on 
y1 = bsplinederiv(4,0,4,x,2);
plot(x,y1);
hold on 
y3 = bsplinederiv(4,0,4,x,3);
plot(x,y3);
hold on 
y4 = bsplinederiv(4,0,4,x,4);
plot(x,y4);
legend('Cubic B-spline','First Derivative','Second Derivative' ,'Third Derivative' ,'Fourth Derivative','location','northwest')
title('First to Fourth Derivatives of Cubic B-spline')