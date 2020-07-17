function plot_surface(Theta,Z,knots,grid_limits,basis_type)
% Plot the surface that comprises the grid of b-splines with corresponding
% coefficients
% Input:
% Theta(Lx1)        - scaling coefficients
% Z(2x2)            - the width of Gaussian function (not important for splines)
% knots(2x2L)       - the support vertexes of each basis function
% grid_limints(4x1) - the region of interes [x1 y1 x2 y2]
% basis_type        - gaussian or bspline
% Distance incriments
dx = (grid_limits(3) - grid_limits(1))/100;%xval
dy = (grid_limits(4) - grid_limits(2))/100;%yval
coordinate_x = [grid_limits(1):dx:grid_limits(3)];
coordinate_y = [grid_limits(2):dy:grid_limits(4)];
N = length(coordinate_x);
M = length(coordinate_y);
[X_grid, Y_grid] = meshgrid(coordinate_x,coordinate_y);
% Suface of basis function superposition
switch basis_type
    case 'gaussian'
        ll = size(knots,2);
        for i = 1:N
            for j = 1:M
                for index1 = 1:ll
                   S = [coordinate_x(i); coordinate_y(j)];
                   c = knots(:,index1);
                   pow = (S - c)'*inv(Z)*((S - c));
                   z(j,i,index1) = Theta(index1)*exp(-pow/2);   
                end
                Z_plot(j,i) = sum(z(j,i,:)) + 1;
            end
        end
    case 'bspline'
         ll = size(knots,2) - 1; %31
%          for t = Theta(1,1:end)
         indexTheta = 1; % index of theta corresponding to each spline
         for index = 1:2:ll %odd numbers
                   support_x = knots(1,index:index+1); %  the x start and ends 
                   support_y = knots(2,index:index+1); %  the y start and ends 
                   Z = tensorproductbspline(4,support_x(1,1),support_x(1,2),support_y(1,1),support_y(1,2),X_grid,Y_grid);
                   Z_plot = Theta(indexTheta)*Z;   
                   Spline(:,:,indexTheta) = Z_plot;
                   indexTheta = indexTheta + 1;
%                    Z_plot = Z*t;
%                   surf(X_grid,Y_grid,Z_plot)
%                   hold on
%          end
         end
       Z_sum = sum(Spline,3);
         
        
     %end
         %for i = 1:N % length of xval (coordinate_x)
           % for j = 1:M % length of yval (coordinate_y)
               % index1 = 1;
               % for index = 1:2:ll %odd numbers
                 %  support_x = knots(1,index:index+1); %  the x start and ends 
                 %  support_y = knots(2,index:index+1); %  the y start and ends 
             
                 %  bf = biorthogonal_spline(coordinate_x(i)/coef_x,coordinate_y(j)/coef_y,support_x/coef_x,support_y/coef_y)
                  % z(j,i,index1) = Theta(index1)*bf 
                 
                 %  index1 = index1 + 1;
               % end
               % Z_plot(j,i) = sum(z(j,i,:)) + 1
        
end
% Z_min = min(Z_sum);
% Z_min_min = min(Z_min);
% if Z_min_min < 0 
%     Z_plot = Z_plot + abs(Z_min_min)*ones(M,N);
% end
% contourf(X_grid,Y_grid,Z_plot,30,'LineWidth',0.5); hold on;
p = surf(X_grid,Y_grid,Z_sum); hold on;
shading interp
alpha(p,0.4);
a_min = min(min(Z_plot));
a_max = max(max(Z_plot));
caxis([a_min a_max]); % the limits of colourbar if it needs to be used. colourbar should be called outside of the function