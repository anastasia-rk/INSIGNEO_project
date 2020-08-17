x = [0:0.2:8];
y = bsplinexval(4,0,2,x);
   p1 =  plot(x,y*2.2,'b');
    hold on 
y1 = bsplinexval(4,1,3,x);
    plot(x,y1*2.23,'b');
y2 = bsplinexval(4,2,4,x);
    plot(x,y2*2,'b');
    y22 = bsplinexval(4,3,5,x);
    plot(x,y22*1.73,'b');
y222 = bsplinexval(4,4,6,x);
    plot(x,y222*1.66,'b');
y3 = bsplinexval(4,5,7,x);
    plot(x,y3*1.87,'b');
y4 = bsplinexval(4,6,8,x);
    plot(x,y4*2.16,'b');
y5 = bsplinexval(4,7,9,x);
    plot(x,y5*2.26,'b');

 y7 = 0.2*sin(x)+1.3
 p2 = plot(x,y7,'r')
 xlabel('x')
 ylabel('f(x)')
 legend([p1 p2],'Fine grid','Coarse grid')
 title('Example of how multigrid captures fast and slow changes in surfaces')

    