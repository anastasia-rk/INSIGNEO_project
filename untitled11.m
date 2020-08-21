folderName='Recruitment/normal_injury/';
iFish = 1;


    load([folderName,'tracks_' num2str(iFish)]);
% Load brightfield image
A = imread([folderName,'bf_' num2str(iFish),'.png']); % for huttenlocher injury 
Cnt = rgb2gray(A); 
[y_max,x_max,zz] = size(A);
% Load the mask
BW = fishmask(Cnt);
% For data split into hours
hour = 0;
if hour~=0
    X = XX;
    Y = YY;
    nTracks = nTracks_t;
end
Tracks = 1:nTracks;
% Create a frame around the image to extend basis function support;
padH = 0; % vertical padding
padW = 200; % horizontal padding
A = padarray(A,[padH, padW]); % creating a padded image
for i=Tracks
    X{i}(:,1) = X{i}(:,1) + padW;
    X{i}(:,2) = X{i}(:,2) + padH;
    Y{i}(:,1) = Y{i}(:,1) + padW;
    Y{i}(:,2) = Y{i}(:,2) + padH;
end
[y_max,x_max,zz] = size(A);
x_lim = [padW x_max-padW];
y_lim = [padH y_max-padH];

   AA = double(BW); % create a surface
Xx = 1:1:size(A,2); % create the grid of x coords
Yy = 1:1:size(A,1); % create the grid of y coords
[Yy_grid,Xx_grid] = meshgrid(Xx,Yy); % mesh
white=[1,1,1]; % surface colour
gg = [0.8,0.8,0.8]; % extra colour for cells 
    

%view(3)
%colormap(my_map);
for j=1:length(Theta)
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
    Theta_temp = zeros(length(Theta),1);                              % zero scaling coeffs for all bfs
    Theta_temp(j) = Theta(j);                                         % only choose one scaling coeff to non-zero
    plot_surface(Theta_temp,Z,knots,grid_limits,basis_type,0.25);
end
hold on
for j=1:length(Theta1)
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
    Theta_temp1 = zeros(length(Theta1),1);                              % zero scaling coeffs for all bfs
    Theta_temp1(j) = Theta1(j);  
    plot_surface(Theta_temp1,Z,knots1,grid_limits,basis_type,1);
end

hold off
