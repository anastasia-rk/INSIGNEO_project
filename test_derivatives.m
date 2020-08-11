clc; clear; local_init;
%% Testing analytical expressions
x = [0:0.1:4];
y = bsplinederiv(4,0,4,x,1);
figure; plot(x,y);

% partial derivative w.r.t to x
for i=1:length(x)
    y(i,:) = [0:0.1:4];
    x_plot(:,i) = x(i);
    dfdx(i,:) = partial_derivative(4,0,4,0,4,x(i),y(i,:),1);
end
figure; surf(x_plot',y',dfdx');
%colormap(my_map);
xlabel('x'); ylabel('y');zlabel('df/dx')
title('Partial Derivative With Respect to X of Cubic B-spline')
clear y x

% partial derivative w.r.t. to y
y = [0:0.1:4];
for i=1:length(y)
    x(i,:) = [0:0.1:4];
    y_plot(:,i) = y(i);
    dfdy(i,:) = partial_derivative(4,0,4,0,4,y(i),x(i,:),1);
end
figure; surf(x,y_plot,dfdy);
%colormap(my_map);
xlabel('x'); ylabel('y');zlabel('df/dy')
title('Partial Derivative With Respect to Y of Cubic B-spline')
