local_init;
%% Select the pattern
nTracks = 100;                                                              % number of tracks to generate in each map
pattern = questdlg('Starting positions distribution', ...
    'Starting Positions',...
	'Uniform','Line','Point','');
switch pattern
    case 'Uniform'
         sx_init = 500*rand(nTracks,1) + 300;
         sy_init = 500*rand(nTracks,1) + 250;
         folderName = 'Simulated/uniform_start/';
    case 'Line'
         sx_init = 450*ones(nTracks,1);
         sy_init = 500*rand(nTracks,1) + 250;
         folderName = 'Simulated/line_start/';
    case 'Point'
         sx_init = 350*ones(nTracks,1);
         sy_init = 350*ones(nTracks,1); 
         folderName = 'Simulated/point_start/';
end
make_folder(folderName);
%% Define field model parameters
basis_type = 'bspline'; % 'gaussian'; %
% Set up limits of the grid: x_min,y_min,x_max,y_max
grid_limits = [0, 0, 1000, 1000];
% Set up number of basis functions
nx = 4; ny = 4; order = 4;
[knots] = setup_spline_support(grid_limits,nx,ny,order);
ll = size(knots,2)/2;
Theta = ones(ll,1);
Theta(1:4,1) = 10;
Theta(5:8,1) = 100;
Theta(9:12,1) = 190;
Theta(13:16,1) = 290;
Theta_model = Theta;
mu = 1;
 %% State-space model parameters
T       = 1;                                                                % sampling time
thta_cv = 0.3;                                                              % reversion to mean in O-U process 
thta_rw = 0.7;                                                              % reversion to mean in O-U process 
I       = eye(2,2);
O       = zeros(2,2);
% mean_vel = 0; % mu of O-U process
% For resonsive mode
F_cv = [I   T*I;...
        O   I- thta_cv*T*I];
% For random walk mode (with I for second order, with O for first order)
F_rw = [I   T*I;...
        O   I - thta_rw*T*I]; 
% Contrl matrix
B_cv = [T^2/2*I; T*I];
B_rw = [O; O];
% Measurement matrix
H = [I O];
% Disturbance matrix
G_cv = [T^2/2*I; T*I];
G_rw = [T^2/2*I; T*I]; 
%% Covariances for the Kalman filter
sig1_Q = 0.5; % RW  disturbance - cell speed
sig2_Q = 0.2; % CV disturbance - random component of cell acceleration
% Process noise - for random walk
Q_rw = sig1_Q*I;  
% Process noise - for responsive mode
Q_cv = sig2_Q*I;
% Measurement noise covariance (mode-independent)
R = 1*eye(2); 
%% Set up Markov chains
n_models = 2;
P_tr = [0.9 0.1;...
        0.1 0.9];...       
F{1} = F_cv;
F{2} = F_rw;
B{1} = B_cv;
B{2} = B_rw;
Q{1} = Q_cv;
Q{2} = Q_rw;
G{1} = G_cv;
G{2} = G_rw;
%%
fig('Heatmap',visFlag); 
colormap(my_map);
Z = 0;
plot_heatmap(Theta_model,Z,knots,grid_limits,basis_type);
colorbar;
hold on;
side = knots(1,2) - knots(1,1);
for i=1:length(Theta_model)
    a = knots(1,i*2-1);
    b = knots(2,i*2-1);
    p=rectangle('Position',[a,b,side,side],'Curvature',0.1,'EdgeColor','w'); 
    hold on;
end
    xlabel('X, a.u.'); ylabel('Y, a.u.');
    print([FigFolder,'Generated_field_',pattern],saveFormat)
%% Modelling cell tracks
for iModel = 1:100
    if disFlag
        iModel
    end
Track_length = 100*ones(nTracks,1); % + floor(100*rand(nTracks,1));
for j=1:nTracks
   % initialise mode
   Mode{j} = zeros(1,Track_length(j));
   r = rand;
   if (r <= .5)
        m = 1;
   else
        m = 2;
   end
   Mode{j}(1) = m;
   chain(1) = m;
   for i=2:Track_length(j)
        this_step_distribution = P_tr(chain(i-1),:);
        cumulative_distribution = cumsum(this_step_distribution);
        r = rand;
        chain(i) = find(cumulative_distribution>r,1); 
        Mode{j}(i) = chain(i);%
   end
   % initialise track    
   X{j} = zeros(4,Track_length(j));
   X{j}(1,1) = sx_init(j,1);
   X{j}(2,1) = sy_init(j,1);
   X{j}(3,1) = 0;
   X{j}(4,1) = 0; 
   % gradient vector
   noise(:,1) = diag(normrnd(0,sqrt(R)));
   grad_model{j} = zeros(Track_length(j),2,ll);
   for i=2:Track_length(j)-1
       m = Mode{j}(i-1);
       % Simulate state
       clear beta;
       grad_model{j}(i-1,:,:) = gradient_bspline(X{j}(1:2,i-1),knots,order);
       aa = F{m}*X{j}(:,i-1);
       beta(:,:) = grad_model{j}(i-1,:,:);
       bb{j}(:,i) = beta*Theta;
       noise_x = mvnrnd([0,0],Q{m})';
       nois_x{j}(:,i) = noise_x;
       X{j}(:,i) = F{m}*X{j}(:,i-1) + B{m}*mu*bb{j}(:,i) + G{m}*noise_x;
       noise(:,i) = mvnrnd([0,0],R)';
       if X{j}(1,i) > 900 || X{j}(1,i) < 0 || X{j}(2,i) > 900 || X{j}(2,i) < 0
           break;
       end
       % Simulate mode
%        probs = P_tr(:,m);
%        r = rand;
%        temp = 0;       
%        for k=1:length(probs)
%            if (r <= probs(k) + temp)
%                m = k;
%                break;
%            else
%                temp = temp +  probs(k);
%            end
%        end
%       Mode{j}(i) = chain(i); 
%        m          = chain(i);
   end
   Xx{j}(:,:) = X{j}(:,1:i);
   clear X{j}
   X{j} = Xx{j};
   Y{j} = H*X{j} + noise(:,:);
   clear noise Xx 
end
%% Example of tracks
if iModel == 1
    fig('Tracks',visFlag); 
    colormap(my_map);
    plot_heatmap(Theta,Z,knots,grid_limits,basis_type);
    hold on;
    for j = 1:nTracks
        txt = num2str(j);
        text(X{j}(1,1)-2,X{j}(2,1), txt,'Color','r','FontSize',15)
        plot(X{j}(1,:),X{j}(2,:),'k','LineWidth',2); hold on; 
        plot(X{j}(1,1),X{j}(2,1),'*r','LineWidth',2); hold on;    
    end
    xlabel('X, a.u.'); ylabel('Y, a.u.');
    print([FigFolder,'Example_tracks_',pattern],saveFormat)
    tikzName = [TikzFolder,'Exampple_trakcs_',pattern,'.tikz'];
    cleanfigure;
    matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '6cm', 'width','6cm','checkForUpdates',false); 
end
%% Save tracks
X_model = X;
Y_model = Y;
clear X Y
for j=1:nTracks
    X{j} = X_model{j};
    Y{j} = Y_model{j}; 
    Mode_model{j} = Mode{j};
end
filename = [folderName,'simulated_tracks_' num2str(iModel)];
save(filename,'X','Y','Mode_model','nTracks','Theta_model');
clear X Y Mode_model 
end