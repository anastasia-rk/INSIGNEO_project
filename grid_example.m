local_init;
%% Setup parameters
basis_type = questdlg('Type of basis function', ...
    'Basis type',...
	'gaussian','bspline','');
% Uncomment the right one:
% basis_type = 'gaussian';
% basis_type = 'bspline';
% Set up limits of the grid: x_min,y_min,x_max,y_max
grid_limits = [0, 0, 1000, 1000];
% Set up number of basis functions along each axis
nx = 4; ny = 4; order = 3;
switch basis_type
    case 'gaussian'                                                         % Symmetric Gaussian function
        [knots,sigma] = setup_gaussian_support(grid_limits,nx,ny);
        % equal sigmas for isotopic basis functions
        Z = [sigma^2 0; 0 sigma^2];                                         % width of the function 
    case 'bspline'
        [knots] = setup_spline_support(grid_limits,nx,ny,order);
        Z = 0;
        % Determine the length of support of a spline
        side1 = knots(1,2) - knots(1,1);                                    % along x axis
        side2 = knots(2,2) - knots(2,1);                                    % along y axis
end
% Arbitrary the scaling coefficients: vector of size (1xL), L = nx*ny
ll = nx*ny;                                                                 % total number of bfs
Theta_model = 10*[1:ll];                                                    % corresponding scaling coeffs
%% Plot the basis funcion grid in 3d
fig('Grid','On'); 
colormap(my_map);
for j=1:length(Theta_model)
   switch basis_type
        case 'gaussian'
            xx = knots(1,j);
            yy = knots(2,j);
        case 'bspline'
            a = knots(1,j*2-1);
            b = knots(2,j*2-1);
            xx = a + side1/2 - 30;
            yy = b + side2/2 - 30;
    end
    Theta_temp = zeros(length(Theta_model),1);                              % zero scaling coeffs for all bfs
    Theta_temp(j) = Theta_model(j);                                         % only choose one scaling coeff to non-zero
    plot_surface(Theta_temp,Z,knots,grid_limits,basis_type);
    hold on;
    txt = (['$\theta_{' num2str(j) '} =$' num2str(Theta_temp(j))]);
    hh = Theta_temp(j);
    text(xx,yy,hh+5,txt,'Color','k','FontSize',18)
end
% view(3)
az = -30;
el = 40;
view(az, el);
zlim([0 max(Theta_model)+20]);
set(gca,'ZTick',[])
set(gca,'XTick',[])
set(gca,'YTick',[])
zlabel('Grid of basis functions')
%% PLot the resultant surface
fig('Surface','On'); 
colormap(my_map);
plot_surface(Theta_model,Z,knots,grid_limits,basis_type);
alpha(1) % overqrite opacity from the function
colorbar;
 %% PLot the gradient
fig('Gradient','On'); 
colormap(my_map);
plot_gradient(Theta_model,Z,knots,grid_limits,basis_type);
%% Plot heatmap
fig('Heatmap','On'); 
colormap(my_map);
plot_heatmap(Theta_model,Z,knots,grid_limits,basis_type);
colorbar; %('north','Color','k','FontSize',12);