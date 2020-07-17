

dx = (grid_limits(3) - grid_limits(1))/100;%xval
dy = (grid_limits(4) - grid_limits(2))/100;%yval
coordinate_x = [grid_limits(1):dx:grid_limits(3)];
coordinate_y = [grid_limits(2):dy:grid_limits(4)];
N = length(coordinate_x);
M = length(coordinate_y);
[X_grid, Y_grid] = meshgrid(coordinate_x,coordinate_y);
ll = size(knots,2) - 1; %31
         %for i = 1:N % length of xval (coordinate_x)
            %for j = 1:M % length of yval (coordinate_y)
                index1 = 1;
                 k = size(Theta_model,2)
                for theta = Theta_model(1,1:k)
                for index = 1:2:ll %odd numbers
                   support_x = knots(1,index:index+1); %  the x start and ends 
                   support_y = knots(2,index:index+1); %  the y start and ends 
                   t = index*theta;
                   Z_plot = tensorproductbspline(3,support_x(1,1),support_x(1,2),support_y(1,1),support_y(1,2),X_grid,Y_grid,t);
                   %coef_x = (support_x(2)-support_x(1))/4;
                   %coef_y = (support_y(2)-support_y(1))/4;
                   %bf = biorthogonal_spline(coordinate_x(i)/coef_x,coordinate_y(j)/coef_y,support_x/coef_x,support_y/coef_y);
                   %z(j,i,index1) = Theta(index1)*bf;  
                   %index1 = index1 + 1;
                   %ZZ = Z_plot*Theta(index1)
                  surf(X_grid,Y_grid,Z_plot)
                  hold on
                end
                end
            
   
         
               