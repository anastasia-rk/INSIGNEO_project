%% Load data
list = {'Normal','Mild','Hutt','Severe'};
[selected,tf] = listdlg('ListString',list,'SelectionMode','single');
Injury = list{selected};
switch Injury
    case 'Normal'
         folderName = 'Recruitment/normal_injury/';
         cc = 1.71; % scaling coefficient for the scale bar
         T  = 2;    % time increment
         nFish = 6; % number of fish to process
    case 'Mild'
         folderName = 'Recruitment/mild_injury/';
         cc = 0.93; % scaling coefficient for the scale bar
         T  = 1.5;  % time increment
         nFish = 5; % number of fish to process
    case 'Hutt'
         folderName = 'Recruitment/huttenlocher_injury/';
         cc = 0.99; % scaling coefficient for the scale bar
         T  = 0.5;  % time increment
         nFish = 3; % number of fish to process
    case 'Severe'
         folderName = 'Recruitment/severe_injury/';
         cc = 1; % scaling coefficient for the scale bar
         T  = 2; % time increment
         nFish = 2; % number of fish to process
end
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
    
figure;
imshow(A); hold on;
% hold on;
for j = Tracks
   plot(Y{j}(:,1),Y{j}(:,2),'-w','LineWidth',1); hold on; 
end
%%%%%%%%%%%%%%%
%for j=Tracks 
%    for i=0:100
%bx=Y{j}(:,1)(abs(Y{j}(:,1)-i)<50); 
%by=Y{j}(:,2)(abs(Y{j}(:,2)-i)<50);
%num_b=numel(bx);
%    end
%end 

surf(Yy_grid,Xx_grid,-AA,'FaceColor',white,'EdgeColor',white);
view(2)
xlim(x_lim);ylim(y_lim);
hold on;
line([250,250+100*cc],[y_max-20,y_max-20],[2,2],'Color','k','LineWidth',5);
txt = ('100 \mu m');
text(250,y_max-70, 2,txt,'Color','k','FontSize',20)
set(gca,'Ydir','reverse')


%% trying to split into tiles 
%figure;

%C = mat2tiles(A,[50,100])
%[m,n] = size(C)
%B=cell2mat(C(randperm(m), randperm(n)));
%imshow(B)

imSz = size(A);
patchSz = [200 350];
xIdxs = [1:patchSz(2):imSz(2) imSz(2)+1];
yIdxs = [1:patchSz(1):imSz(1) imSz(1)+1];
patches = cell(length(yIdxs)-1,length(xIdxs)-1)
for i = 1:length(yIdxs)-1
    Isub = A(yIdxs(i):yIdxs(i+1)-1,:);
    for j = 1:length(xIdxs)-1
        patches{i,j} = Isub(:,xIdxs(j):xIdxs(j+1)-1);
    end
end
figure, imagesc(patches{3,3})
hold on
for j = Tracks
   plot(Y{j}(:,1),Y{j}(:,2),'-w','LineWidth',1); hold on; 
end
%% histogram plot for concentrations - didnt really work 
figure;
for j = Tracks
    nbins=[90 90];
[N,C]=hist3([Y{j}(:,1), Y{j}(:,2)],nbins)
contourf(C{1},C{2},N)
  hold on 
end
%%

figure;
basis_type = 'bspline';
% Set up limits of the grid: x_min,y_min,x_max,y_max
grid_limits = [0, 0, x_max, y_max];
% Set up number of basis functions
 nx=8; ny=8; order=4;
[knots] = setup_spline_support(grid_limits,nx,ny,order); % spline support nodes
Z = 0;
ll = size(knots,2)/2; % size of parameter vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
side1 = knots(1,2) - knots(1,1);                             % along x axis
side2 = knots(2,2) - knots(2,1);                             % along y axis
mx=4; my=4
[knots1] = setup_spline_support(grid_limits,mx,my,order);
side11 = knots1(1,2)-knots(1,1);
side22 = knots1(2,2)-knots(2,1);
tt = size(knots1,2)/2;
Theta1 = 1:mx*my; %zeros(tt,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise field model parameters
grid_limits1 = grid_limits;
%Theta = zeros(ll,1);
Theta = 1:nx*ny;
% Propotional coefficient (sensitivity)
mu_field = 1;


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