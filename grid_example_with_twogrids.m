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
nx =  6;ny = 6; order = 4;
mx = 3; my = 3;
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
        
        [knots1] = setup_spline_support(grid_limits,mx,my,order);
        side11 = knots1(1,2)-knots(1,1);
        side22 = knots1(2,2)-knots(2,1);
end
% Arbitrary the scaling coefficients: vector of size (1xL), L = nx*ny
ll = nx*ny;   
tt = mx*my;
Theta_model = 1*ones*[1:ll];    
Theta_model1 = 5*ones*[1:tt]; % ones(ll,1); %                          % corresponding scaling coeffs
% total number of bfs
%% Initiate a field as a flat surface
grid_limits2 = [100, 200, 720, 750];
[Theta] = initiate_field(0,100,grid_limits2,knots);
hold on;
plot_surface(Theta,Z,knots,grid_limits,basis_type,0.5);
%Theta_model = Theta;
%% Plot the basis funcion grid in 3d
fig('Grid','On'); 
ax1 = subplot(2,1,1);
%colormap(my_map);
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
    plot_surface_old(Theta_temp,Z,knots,grid_limits,basis_type,0.25);
    %txt = (['$\theta_{' num2str(j) '} =$' num2str(Theta_temp(j))]);
    %hh = Theta_temp(j);
    %text(xx,yy,hh,txt,'Color','k','FontSize',18)
end
colormap(ax1,spring)
ax2 =  subplot(2,1,2);
%hold on 
for j=1:length(Theta_model1)
   switch basis_type
        case 'gaussian'
            xx = knots1(1,j);
            yy = knots1(2,j);
        case 'bspline'
            a = knots1(1,j*2-1);
            b = knots1(2,j*2-1);
            xx = a + side1/2 - 30;
            yy = b + side2/2 - 30;
            
   end
    Theta_temp1 = zeros(length(Theta_model1),1);                              % zero scaling coeffs for all bfs
    Theta_temp1(j) = Theta_model1(j);  
    plot_surface_old(Theta_temp1,Z,knots1,grid_limits,basis_type,0.5);
    %txt = (['$\theta_{' num2str(j) '} =$' num2str(Theta_temp(j))]);
    %hh = Theta_temp(j);
    %text(xx,yy,hh,txt,'Color','k','FontSize',18)
end
  colormap(ax2,winter)
% view(3)
az = -30;
el = 40;
view(az, el);
zlim([0 max(Theta_model)+20]);
% set(gca,'ZTick',[])
set(gca,'XTick',[])
set(gca,'YTick',[])
zlabel('Grid of basis functions')
%% PLot the resultant surface
fig('Surface','On'); 
%colormap(my_map);
plot_surface(Theta_model,Z,knots,grid_limits,basis_type,1);
% alpha(1) % overqrite opacity from the function
colorbar;
 %% PLot the gradient
fig('Gradient','On'); 
%colormap(my_map);
plot_gradient(Theta_model,Z,knots,grid_limits,basis_type);
%% Plot heatmap
fig('Heatmap','On'); 
%colormap(my_map);
plot_heatmap(Theta_model,Z,knots,grid_limits,basis_type);
colorbar; %('north','Color','k','FontSize',12);