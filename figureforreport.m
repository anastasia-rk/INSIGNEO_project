x = [0:0.2:12];
y = bsplinexval(4,0,2,x);
   p1 =  plot(x,y*2.2,'b');
    hold on 
y1 = bsplinexval(4,2,4,x);
    plot(x,y1*2,'b');
y2 = bsplinexval(4,4,6,x);
    plot(x,y2*1.65,'b');
y3 = bsplinexval(4,6,8,x);
    plot(x,y3*2.15,'b');
y4 = bsplinexval(4,8,10,x);
    plot(x,y4*2.08,'b');
y5 = bsplinexval(4,10,12,x);
    plot(x,y5*1.66,'b');

 y7 = 0.2*sin(x)+1.3
 p2 = plot(x,y7,'r')
 xlabel('x')
 ylabel('f(x)')
 legend([p1 p2],'Fine grid','Coarse grid')
 title('Example of how multigrid captures fast and slow changes in surfaces')

 A = [1 2 3;4 5 6]
 reshape(A,1,[])'
    