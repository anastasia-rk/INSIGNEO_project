function plot_heatmap(Theta,Z,knots,grid_limits,basis_type)
% Plot the surface that comprises the grid of b-splines with corresponding
% coefficients
% Input:
% Theta(Lx1)        - scaling coefficients
% Z(2x2)            - the width of Gaussian function (not important for splines)
% knots(2x2L)       - the support vertexes of each basis function
% grid_limints(4x1) - the region of interes [x1 y1 x2 y2]
% basis_type        - gaussian or bspline
% Distance incriments
dx = (grid_limits(3) - grid_limits(1))/100;
dy = (grid_limits(4) - grid_limits(2))/100;
coordinate_x = [grid_limits(1):dx:grid_limits(3)];
coordinate_y = [grid_limits(2):dy:grid_limits(4)];
N = length(coordinate_x);
M = length(coordinate_y);
[X_grid, Y_grid] = meshgrid(coordinate_x,coordinate_y);
% Suface of spline superposition
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
                Z_plot(j,i) = sum(z(j,i,:));
            end
        end
    case 'bspline'
         ll = size(knots,2) - 1;
         indexTheta = 1; % index of theta corresponding to each spline
         for index = 1:2:ll %odd numbers
                   support_x = knots(1,index:index+1); %  the x start and ends 
                   support_y = knots(2,index:index+1); %  the y start and ends 
                   Z = tensorproductbspline(4,support_x(1,1),support_x(1,2),support_y(1,1),support_y(1,2),X_grid,Y_grid);
                   Z_plot = Theta(indexTheta)*Z;   
                   Spline(:,:,indexTheta) = Z_plot;
                   indexTheta = indexTheta + 1;                 
         end
         Z_sum = sum(Spline,3);
end
Z_min = min(Z_sum);
Z_min_min = min(Z_min);
a_min = min(min(Z_sum-Z_min));
a_max = max(max(Z_sum-Z_min));
p = pcolor(X_grid,Y_grid,Z_sum-Z_min);
view(2);
shading interp
caxis manual
alpha(p,0.8);
caxis([a_min a_max]);