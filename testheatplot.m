folderName='Recruitment/normal_injury/';
iFish = 2;
% Load tracking data
load([folderName,'tracks_' num2str(iFish)]);
% Load brightfield image
A = imread([folderName,'bf_' num2str(iFish),'.png']); % for huttenlocher injury 
Cnt = rgb2gray(A); 
[y_max,x_max,zz] = size(A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BW = fishmask(Cnt);
%%
padH = 000; % vertical padding
padW = 200;% horizontal padding
A = padarray(A,[padH, padW]); % creating a padded image
for i=1:nTracks
    X{i}(:,1) = X{i}(:,1) + padW;
    X{i}(:,2) = X{i}(:,2) + padH;
    Y{i}(:,1) = Y{i}(:,1) + padW;
    Y{i}(:,2) = Y{i}(:,2) + padH;
end
[y_max,x_max,zz] = size(A);
x_lim = [padW x_max-padW];
y_lim = [padH y_max-padH];
%
%%
AA = double(BW); % create a surface
Xx = 1:1:size(A,2); % create the grid of x coords
Yy = 1:1:size(A,1); % create the grid of y coords
[Yy_grid,Xx_grid] = meshgrid(Xx,Yy); % mesh
white=[1,1,1]; % surface colour
gg = [0.8,0.8,0.8]; % extra colour for cells
figure; set(gcf,'color','w');
imshow(A); hold on;
% hold on;
% %% Resample tracks according to the most proable mode
% for k=Tracks
%     for t=1:T(k)
%         P_for_resampling{k}(t,:,:) = P_out{k}{t};
%     end
% end
% nTracks = length(Tracks);
% [X_resample,P_resample,nResampled,ind_resample] = resample_tracks(X_out,P_for_resampling,mode_probable,1,nTracks);
% [X_rw,P_rw,nRW,ind_rw] = resample_tracks(X_out,P_for_resampling,mode_probable,1,nTracks);
% [X_dead,P_dead,nDead,ind_dead] = resample_tracks(X_out,P_for_resampling,mode_probable,1,nTracks);
% for j = Tracks
%    plot(X_out{j}(1,:),X_out{j}(2,:),'-b','LineWidth',1); hold on; 
%    tx = num2str(j);
%    text(X_out{j}(1,1),X_out{j}(2,1),2, tx, 'interpreter', 'latex','Color','k','FontSize',12);
% end
%  for j = 1:nResampled
%    plot(X_resample{j}(1:end-1,1),X_resample{j}(1:end-1,2),'-g','LineWidth',1); hold on; 
%  end
%   for j = 1:nDead
%    plot(X_dead{j}(1:end-1,1),X_dead{j}(1:end-1,2),'-','Color',white,'LineWidth',1); hold on; 
%   end
% %  for j = Tracks % mark the starting point of each track
% %    if Mode{j}(1) == 1
% %         plot(X_out{j}(1,1),X_out{j}(2,1),'*g','LineWidth',1); hold on;
% %    elseif Mode{j}(1) == 2
% %         plot(X_out{j}(1,1),X_out{j}(2,1),'*b','LineWidth',1); hold on;
% %    else
% %         plot(X_out{j}(1,1),X_out{j}(2,1),'*','Color',white,'LineWidth',1); hold on;
% %    end
% %  end
% xlim(x_lim);ylim(y_lim);
% hold on;
% surf(Yy_grid,Xx_grid,-AA,'FaceColor',white,'EdgeColor',white);
% view(2)
% hold on;
% hold on;
% text(625+ padW, 825 + padH,2, '$\star$  ', 'interpreter', 'latex','Color','g','FontSize',20);
% text(650+ padW, 825 + padH,2, ' - $M^1$', 'interpreter', 'latex','Color','k','FontSize',20);
% text(750+ padW, 825 + padH,2, '$\star$  ', 'interpreter', 'latex','Color','b','FontSize',20);
% text(775+ padW, 825 + padH,2, ' - $M^2$', 'interpreter', 'latex','Color','k','FontSize',20);
% text(875+ padW, 825 + padH,2, '$\star$  ', 'interpreter', 'latex','Color',white,'FontSize',20);
% text(900+ padW, 825 + padH,2, ' - $M^3$', 'interpreter', 'latex','Color','k','FontSize',20);
% set(gca,'Ydir','reverse')
% % line([50,221],[825,825],'Color','k','LineWidth',5);
%  line([250,250+100*cc],[y_max-20,y_max-20],[2,2],'Color','k','LineWidth',5);
% % line([50,221],[575,575],'Color','k','LineWidth',5);
% txt = ('100 $\mu$m');
% text(250,y_max-45, 2,txt,'Color','k','FontSize',20)
% %text(70,550,2,txt,'Color','k','FontSize',20)
% % set(gca,'Ydir','reverse')
%% Field model parameters
%basis_type = 'gaussian';
basis_type = 'bspline';
% Set up limits of the grid: x_min,y_min,x_max,y_max
grid_limits = [0, 0, x_max, y_max];
% Set up number of basis functions
nx = 6; ny = 6; order = 4;
[knots] = setup_spline_support(grid_limits,nx,ny,order); % spline support nodes
Z = 0;
ll = size(knots,2)/2; % size of parameter vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
side1 = knots(1,2) - knots(1,1);                                    % along x axis
side2 = knots(2,2) - knots(2,1);                                    % along y axis

mx = 4; my = 4;
[knots1] = setup_spline_support(grid_limits,mx,my,order);
side11 = knots1(1,2)-knots(1,1);
side22 = knots1(2,2)-knots(2,1);
tt = size(knots1,2)/2;
Theta1 = 1:16 %zeros(tt,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise field model parameters
grid_limits1 = grid_limits;
%Theta = zeros(ll,1);
Theta = 1:36;
% Propotional coefficient (sensitivity)
mu_field = 1;
        
        
        
%fig('Heatmap',visFlag);
imshow(A); hold on;
plot_heatmap(Theta,Z,knots,grid_limits,basis_type);
hold on
plot_heatmap(Theta1,Z,knots1,grid_limits,basis_type);
% alpha(0.5)
hold on;
surf(Yy_grid,Xx_grid,-AA,'FaceColor',white,'EdgeColor',white);
view(2)
hold on;
xlim(x_lim);ylim(y_lim);
line([250,250+100*cc],[y_max-20,y_max-20],[2,2],'Color','k','LineWidth',5);
% line([50,221],[575,575],'Color','k','LineWidth',5);
txt = ('100 $\mu$m');
text(250,y_max-45, 2,txt,'Color','k','FontSize',20)
% text(70,550, 2,txt,'Color','k','FontSize',20)
set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(gca,'Ydir','reverse');
%print([FigFolder,'heatmap_',Injury,num2str(iFish)],saveFormat)

%fig('Surface',visFlag);
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


