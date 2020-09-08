% A.R. Kadochnikova - 14/06/2016

% Simultaneous cell mode identification and chemotactic field inferrence.
% State estimation/mode identification carried out through IMM
% forward-backward smoother
% Filed inferrence carried out through EM algorithm.

% [1] V. Kadirkamanathan, G.R. Holmes "The Neutrophils Eye-View: Inference
% and Visualisation of the Chemoattractant field Driving Cell Chemotaxis", 2012. 
% [2] V. Kadirkamanathan, S. Anderson "Maximum-Likelihood Estimation of
% Delta-Domain Model Parameters from Noisy Output Signals", 2008.
% [3] Kalman Smoother without inv(P) A.Logothesis and V.Krishnamurthy "EM 
% alforithms for MAP estimation of jump Markov linear systems", 1999.
% [4] Multiple model Kalman filter from X. Zhang "Modelling and
% Identification of Neutrophil Cell Dynamic Behaviour", 2016.
local_init;
visFlag = 'Off';
%% Load data
% list = {'Normal','Mild','Hutt','Severe'};
% [selected,tf] = listdlg('ListString',list,'SelectionMode','single');
Injury = 'Normal';
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
         nFish = 4; % number of fish to process
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
%%
iFish = 3;
% Load tracking data
load([folderName,'tracks_' num2str(iFish)]);
% Load brightfield image
A = imread([folderName,'bf_' num2str(iFish),'.png']); % for huttenlocher injury 
Cnt = rgb2gray(A); 
[y_max,x_max,zz] = size(A);
% Load the mask
load([folderName,'mask_',num2str(iFish)]); % downloads variable BW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Room for Kate
% use your edge detection function to produce a mask from image A
% The mask variable should be named BW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For data split into hours
hour = 0;
if hour~=0
    X = XX;
    Y = YY;
    nTracks = nTracks_t;
end
%% Create a frame around the image to extend basis function support;
padH = 0; % vertical padding
padW = 200; % horizontal padding
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Room for Kate
% plot the padded image A to see what has changed from the original imgage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Field model parameters
%basis_type = 'gaussian';
basis_type = 'bspline';
% Set up limits of the grid: x_min,y_min,x_max,y_max
grid_limits = [0, 0, x_max, y_max];
% Set up number of basis functions
nx = 5; ny = 4; order = 4;
[knots] = setup_spline_support(grid_limits,nx,ny,order); % spline support nodes
Z = 0;
ll = size(knots,2)/2; % size of parameter vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Room for Kate
% create a multi-resolution grid of b-splines
% coars and fine grid.  knots - support points, ll - number of splines (overall)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise field model parameters
grid_limits1 = grid_limits;
Theta = zeros(ll,1);
% Propotional coefficient (sensitivity)
mu_field = 1;
%% State-space model parameters
% Three candidate models of cell dynamics are examined: 
fprintf(1,'Setting model parameters... \n')
x_len = 4; % size of the state vector
% Transition matrix
I =   eye(2,2);
O = zeros(2,2);
thta1 = 0.3; % reversion to mean in O-U process 
thta2 = 0.3; % reversion to mean in O-U process 
thta3 = 0.5; % reversion to mean in O-U process 

mean_vel = 0; % mu of O-U process
% brownian motoion with friction  (velocity as O-U process)
% For CV
F_cv = [I T*I;...
    O   I - thta1*T*I];  
% For RW
F_rw = [I T*I;...
    O   I - thta2*T*I]; 
% For stationary
F_st = [I T*I;...
    O   I - thta2*T*I];  
% Contrl matrix
B_cv = [T^2/2*I; T*I];
B_rw = [O; O];
% Measurement matrix
C = [I O];
% Disturbance matrix
G_cv = [T^2/2*I; T*I];
G_rw = [T^2/2*I; T*I]; 
clear T
%% Disturbance matrices
% Q - describes power of state noise
sig1_Q = 2; % RW  disturbance - cell speed
sig2_Q = 2; % CV disturbance - random component of cell acceleration
sig3_Q = 0.5;
% For DIFF
m = sig1_Q*G_cv*G_cv';  
d = diag(m);
Q_rw = diag(d);
Q_rw_inv = (1/sig1_Q)*eye(2);
% For DRIFT
m = sig2_Q*G_cv*G_cv';
d = diag(m);
Q_cv = diag(d);
Q_cv_inv = (1/sig2_Q)*eye(2);
% For STAT
m = sig3_Q*G_cv*G_cv';
d = diag(m);
Q_st = diag(d);
Q_st_inv = (1/sig3_Q)*eye(2);
% R - measurement noise
sig2_R = 2;
R = sig2_R*eye(2);
clear d m
%% Set up immf
n_models = 3;
p_tr = [0.8 0.1 0.1;...
        0.1 0.8 0.1;...
        0.1 0.1 0.8];  
Candidate = [1:1:n_models];
% 1 - drifting (keller-segel)
% 2 - diffusing (random walk)
% 3 - stationary (random walk)
F{1} = F_cv;
F{2} = F_rw;
F{3} = F_st;
B{1} = B_cv;
B{2} = B_rw;
B{3} = B_rw;
Q{1} = Q_cv;
Q{2} = Q_rw;
Q{3} = Q_st;
G{1} = G_cv;
G{2} = G_rw;
G{3} = G_rw;

Sig_w{1} = G{1}*(inv(G{1}'*G{1}))'*Q_cv_inv*inv(G{1}'*G{1})*G{1}';
Sig_w{2} = G{2}*(inv(G{2}'*G{2}))'*Q_rw_inv*inv(G{2}'*G{2})*G{2}';
Sig_w{3} = G{2}*(inv(G{3}'*G{3}))'*Q_st_inv*inv(G{3}'*G{3})*G{3}';

mu_0 = ones(1,n_models);
mu_0 = mu_0./n_models;

switch_from = zeros(3,1);
switch_to = zeros(3,1);
switches = zeros(3,3);
nSwitch = 0;
%% Estimation framework
Tracks = [1:nTracks]; 
for k=Tracks % Assign track lengths outside of parfor loop for use in sliced variables
    Y{k} = Y{k}';
    T(k) = length(Y{k});
end
% Create parallel pool
converged  = false;
iter_max   = 20;
iter       = 0;
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local');
end
%% EM loop
while  (iter < iter_max)
%% Expectation step
fprintf('E-step... \n')
iter = iter + 1;  
%% Parfor loop for IMM
tic
parfor k=Tracks     
%% recurtion cycle for IMM filter
    [x_m,P,mu,mode_probable_f{k}] = IMM_forward(T(k),x_len,n_models,p_tr,mu_0,Y{k},Q,R,F,B,C,Theta,knots,order);
%% Recurtion cycle for IMM RTS smoother
    [x_s,P_s,mu_s,x_merged,P_merged,mu_cond{k},mode_probable{k}] = IMM_backward(T(k),x_len,n_models,p_tr,mu,x_m,P,F,B,Q,Theta,knots,order);
%% Output of the IMM algorithm    
    [X_out{k},P_out{k},Mu_out{k},X_cond{k},P_cond{k}] = IMM_output(T(k),n_models,x_merged,P_merged,x_s,P_s,mu_s);
end % for track (k)
toc
%% Maximization step [2]
fprintf('M-step... \n')
%% Transition probability matrix 
mu_joint = zeros(n_models,n_models);
for j=1:n_models
    denom(j) = 0;
    for l=1:n_models
        for k=Tracks
            for t=1:T(k)-2
                mu_joint(j,l) = mu_joint(j,l) + Mu_out{k}(j,t)*mu_cond{k}{t+1}(j,l);
            end % for time (t)
        end % for trakcs (k)
        denom(j) = denom(j) + mu_joint(j,l);
    end% for modes (l)
    for l = 1:n_models
        transitions{iter}(j,l) = mu_joint(j,l)/denom(j);
    end % for modes (l)
end % for modes (j)
p_tr_old = p_tr;
p_tr = transitions{iter}
%% Field parameters
sum1 = 0; sum2 =  0; sum3 = 0;
tic
parfor k=Tracks
     for t=1:T(k)-1
         x_tilde = X_out{k}(:,t);
         grad = gradient_bspline(x_tilde(1:2),knots,order);
         for j=1:n_models
             x_j_tilde = X_cond{k}{j}(:,t+1);
             bb = B{j}*mu_field*grad;
             sum1 = sum1 + Mu_out{k}(j,t)*bb'*Sig_w{j}*bb;
             sum2 = sum2 + Mu_out{k}(j,t)*bb'*Sig_w{j}*(x_j_tilde - F{j}*x_tilde);
         end % for modes (j) 
     end % for times (t)
end % for tracks (k)
toc   
% Parameter MLE
Theta_old = Theta;
Theta = pinv(sum1)*sum2;
FieldCoeff{iter} = Theta;
% Expected Fisher information
Expected_info = sum1;
Missing_info = (sum1*Theta - sum2)*(sum1*Theta_old - sum2)';
% temp(:,iter) = Theta;
%% Check convergence
fprintf('Checking convergence condition... \n')
[converged_theta]   = converge_param(Theta,Theta_old,iter);
[converged_phi]     = converge_param(p_tr(:),p_tr_old(:),iter);
% [converged_r]       = converge_param(vec(R),vec(R_old),iter);
    if converged_theta & converged_phi % & converged_r                      % if convergence condition fulfilled
        break;
    end
end % for EM iterations (iter)
%% Resample tracks according to the most proable mode
for k=Tracks
    for t=1:T(k)
        P_for_resampling{k}(t,:,:) = P_out{k}{t};
    end
end
nTracks = length(Tracks);
[X_resample,P_resample,nResampled,ind_resample] = resample_tracks(X_out,P_for_resampling,mode_probable,1,nTracks);
[X_rw,P_rw,nRW,ind_rw] = resample_tracks(X_out,P_for_resampling,mode_probable,1,nTracks);
[X_dead,P_dead,nDead,ind_dead] = resample_tracks(X_out,P_for_resampling,mode_probable,1,nTracks);
%% Save data to file
if hour==0
  filename = ['Results/',Injury,num2str(iFish)];
else
  filename = ['Results/',Injury,num2str(iFish),'hr_',num2str(hour)]; 
end
save(filename);
%% Confidence region volume
Observed_info = Expected_info - Missing_info;
COV_ALL = pinv(Observed_info);
VOL_ALL = ellipsoid_volume(ll,0.05,COV_ALL);
VOL_ALL
TOTALV = trace(COV_ALL)/length(Theta);
TOTALV
%% Count mode transistions - see how it compares with the estimate p_tr
for k=Tracks
    [nn,T] = size(mode_probable{k});  
 for t = 2:T
        switch_to(mode_probable{k}(t)) = switch_to(mode_probable{k}(t)) + 1;
        switches(mode_probable{k}(t-1),mode_probable{k}(t)) = switches(mode_probable{k}(t-1),mode_probable{k}(t)) + 1;
 end
end
nSwitch = nSwitch +  sum(switch_to);

prob_freq = zeros(3);
for i=1:3
   prob_freq(:,i) =  switches(:,i)/switch_to(i);
end
%% Transition probability matrix for all models - statistics
iter_plot = iter+1;
for j=1:n_models
    for l=j:n_models
        phi_mean{j,l}(1) = 0.5;
        for it = 2:iter_plot
            phi_mean{j,l}(it) = mean(transitions{it-1}(j,l));
        end
    end
end
fig('Convergence of Phi',visFlag); 
% line([0, iter_plot],[0.9, 0.9],'Color','black'); hold on;
% line([0, iter_plot],[0.8, 0.8],'Color','black'); hold on;
% line([0, iter_plot],[0.2, 0.2],'Color','black'); hold on;
% line([0, iter_plot],[0.1, 0.1],'Color','black'); hold on;
x_patch = [0:iter];
iPhi = 0;
for j=1:n_models
    for l=j:n_models
        iPhi = iPhi + 1;
        thisline(iPhi) = plot(x_patch,phi_mean{j,l}); hold on;
        names{iPhi} = ['$\hat{\phi}_{',num2str(j),num2str(l),'}$'];
%         color = get(thisline(iPhi), 'Color');
%         yy = [phi_mean{j,l}+3*phi_std{j,l} fliplr(phi_mean{j,l}-3*phi_std{j,l})];
%         patch(xxx,yy,color,'EdgeColor',color,'FaceAlpha',0.2); hold on;
%         for it=1:iter_plot
%             scatter(it*ones(1,iModel),phi_est{j,l}(it,:),'filled','MarkerFaceColor',color); hold on %
%         end
%         clear color
    end
end
ylim([0,1]);
xlabel('Iteration of the EM algorithm');
ylabel('$\hat{\phi}_{j,i}$')
legend(thisline,names);
% {'$\hat{\phi}_{1,1}$','$\hat{\phi}_{1,2}$','$\hat{\phi}_{2,1}$','$\hat{\phi}_{2,2}$'}
print([FigFolder,'phi_convergence_',Injury,num2str(iFish)],saveFormat)
matlab2tikz([TikzFolder,'phi_convergence_',Injury,num2str(iFish),'.tikz'], 'showInfo', false,'parseStrings',false, ...
         'standalone', false,'height', '3cm', 'width','4cm');
%% Extract velocities from resampled tracks
Vx_count = 1;
Vy_count = 1;
for j = Tracks
   [n,l] = size(X_out{j});
   Vx(1,Vx_count:Vx_count+l-1) = X_out{j}(3,:)/cc;
   Vy(1,Vy_count:Vy_count+l-1) = X_out{j}(4,:)/cc;
   Vx_count = Vx_count + l;
   Vy_count = Vy_count + l;
end
Vx_count = 1;
Vy_count = 1;
for j = 1:nResampled
   [l,n] = size(X_resample{j});
   Vx_1(1,Vx_count:Vx_count+l-1) = X_resample{j}(:,3)'./cc;
   Vy_1(1,Vy_count:Vy_count+l-1) = X_resample{j}(:,4)'./cc;
   Vx_count = Vx_count + l;
   Vy_count = Vy_count + l;
end
Vx_count = 1;
Vy_count = 1;
for j = 1:nRW
   [l,n] = size(X_rw{j});
   Vx_2(1,Vx_count:Vx_count+l-1) = X_rw{j}(:,3)'./cc;
   Vy_2(1,Vy_count:Vy_count+l-1) = X_rw{j}(:,4)'./cc;
   Vx_count = Vx_count + l;
   Vy_count = Vy_count + l;
end
Vx_count = 1;
Vy_count = 1;
for j = 1:nDead
   [l,n] = size(X_dead{j});
   Vx_3(1,Vx_count:Vx_count+l-1) = X_dead{j}(:,3)'./cc;
   Vy_3(1,Vy_count:Vy_count+l-1) = X_dead{j}(:,4)'./cc;
   Vx_count = Vx_count + l;
   Vy_count = Vy_count + l;
end
%% Velocity histogramms
numberOfBins = 50;
Med_x = median(Vx_1);
Mean_x = mean(Vx_1);
Med_y = median(Vy_1);
Mean_y = mean(Vy_1);
yylim = [0 55];
fig('Vx hist. 1',visFlag)
% title('Velocity histogram')
[counts, binValues] = hist(Vx_1, numberOfBins);
normalizedCounts = 100 * counts / sum(counts);
p1 = bar(binValues, normalizedCounts, 'barwidth', 1);
set(p1,'FaceColor','white');
hold on;
line([Med_x, Med_x], yylim,'LineStyle','-.', 'LineWidth', 1, 'Color', 'k');
line([Mean_x, Mean_x], yylim, 'LineWidth', 1, 'Color', 'r');
text(Mean_x +2, 50,['$\mu_x$=',num2str(Mean_x)],'Color','k','FontSize',20);
xlabel('$v_x$, $\mu$m/min', 'interpreter', 'latex');
ylabel('\textrm{$\%$}', 'interpreter', 'latex');
% 
matlab2tikz([TikzFolder,'vx_',Injury,num2str(iFish),'.tikz'], 'showInfo', false,'parseStrings',false, ...
         'standalone', false,'height', '3cm', 'width','4cm');
fig('Vy hist. 1',visFlag)
[counts, binValues] = hist(Vy_1, numberOfBins);
normalizedCounts = 100 * counts / sum(counts);
p1 = bar(binValues, normalizedCounts, 'barwidth', 1);
set(p1,'FaceColor','white');
hold on;
line([Med_y, Med_y], yylim,'LineStyle','-.','LineWidth', 1, 'Color', 'k');
line([Mean_y, Mean_y], yylim, 'LineWidth', 1, 'Color', 'r');
text(Mean_y +2, 50,['$\mu_y$=',num2str(Mean_y)],'Color','k','FontSize',20);
xlabel('$v_y$, $\mu$m/min', 'interpreter', 'latex');
ylabel('\textrm{$\%$}', 'interpreter', 'latex');
matlab2tikz([TikzFolder,'vy1_',Injury,num2str(iFish),'.tikz'], 'showInfo', false,'parseStrings',false, ...
         'standalone', false,'height', '3cm', 'width','4cm');
Med_x = median(Vx_2);
Mean_x = mean(Vx_2);
Med_y = median(Vy_2);
Mean_y = mean(Vy_2);
yylim = [0 55];
fig('Vx hist. 2',visFlag)
% title('Velocity histogram')
[counts, binValues] = hist(Vx_2, numberOfBins);
normalizedCounts = 100 * counts / sum(counts);
p1 = bar(binValues, normalizedCounts, 'barwidth', 1);
set(p1,'FaceColor','white');
hold on;
line([Med_x, Med_x], yylim,'LineStyle','-.', 'LineWidth', 1, 'Color', 'k');
line([Mean_x, Mean_x], yylim, 'LineWidth', 1, 'Color', 'r');
text(Mean_x +2, 50,['$\mu_x$=',num2str(Mean_x)],'Color','k','FontSize',20);
xlabel('$v_x$, $\mu$m/min', 'interpreter', 'latex');
ylabel('\textrm{$\%$}', 'interpreter', 'latex');
matlab2tikz([TikzFolder,'vx2_',Injury,num2str(iFish),'.tikz'], 'showInfo', false,'parseStrings',false, ...
         'standalone', false,'height', '3cm', 'width','4cm');
fig('Vy hist. 2',visFlag)
[counts, binValues] = hist(Vy_2, numberOfBins);
normalizedCounts = 100 * counts / sum(counts);
p1 = bar(binValues, normalizedCounts, 'barwidth', 1);
set(p1,'FaceColor','white');
hold on;
line([Med_y, Med_y], yylim,'LineStyle','-.','LineWidth', 1, 'Color', 'k');
line([Mean_y, Mean_y], yylim, 'LineWidth', 1, 'Color', 'r');
text(Mean_y +2, 50,['$\mu_y$=',num2str(Mean_y)],'Color','k','FontSize',20);
xlabel('$v_y$, $\mu$m/min', 'interpreter', 'latex');
ylabel('\textrm{$\%$}', 'interpreter', 'latex');
matlab2tikz([TikzFolder,'vy2_',Injury,num2str(iFish),'.tikz'], 'showInfo', false,'parseStrings',false, ...
         'standalone', false,'height', '3cm', 'width','4cm');
      
Med_x = median(Vx_3);
Mean_x = mean(Vx_3);
Med_y = median(Vy_3);
Mean_y = mean(Vy_3);
yylim = [0 55];
fig('Vx hist. 3',visFlag)
[counts, binValues] = hist(Vx_3, numberOfBins);
normalizedCounts = 100 * counts / sum(counts);
p1 = bar(binValues, normalizedCounts, 'barwidth', 1);
set(p1,'FaceColor','white');
hold on;
line([Med_x, Med_x], yylim,'LineStyle','-.', 'LineWidth', 1, 'Color', 'k');
line([Mean_x, Mean_x], yylim, 'LineWidth', 1, 'Color', 'r');
text(Mean_x +2, 50,['$\mu_x$=',num2str(Mean_x)],'Color','k','FontSize',20);
xlabel('$v_x$, $\mu$m/min', 'interpreter', 'latex');
ylabel('\textrm{$\%$}', 'interpreter', 'latex');
matlab2tikz([TikzFolder,'vx3_',Injury,num2str(iFish),'.tikz'], 'showInfo', false,'parseStrings',false, ...
         'standalone', false,'height', '3cm', 'width','4cm');
fig('Vy hist. 3',visFlag)
[counts, binValues] = hist(Vy_3, numberOfBins);
normalizedCounts = 100 * counts / sum(counts);
p1 = bar(binValues, normalizedCounts, 'barwidth', 1);
set(p1,'FaceColor','white');
hold on;
line([Med_y, Med_y], yylim,'LineStyle','-.','LineWidth', 1, 'Color', 'k');
line([Mean_y, Mean_y], yylim, 'LineWidth', 1, 'Color', 'r');
text(Mean_y +2, 50,['$\mu_y$=',num2str(Mean_y)],'Color','k','FontSize',20);
xlabel('$v_y$, $\mu$m/min', 'interpreter', 'latex');
ylabel('\textrm{$\%$}', 'interpreter', 'latex');
matlab2tikz([TikzFolder,'vy3_',Injury,num2str(iFish),'.tikz'], 'showInfo', false,'parseStrings',false, ...
         'standalone', false,'height', '3cm', 'width','4cm');

%%
fig('Positions',visFlag);
subplot(2,2,1);
for k=Tracks
    plot(Y{k}(1,:),Y{k}(2,:),'k'); hold on;
    plot(X_cond{k}{1}(1,2:end),X_cond{k}{1}(2,2:end),'r'); hold on;
end
xlabel('X');ylabel('y')
title('Smoothed mode 1')
subplot(2,2,2);
for k=Tracks
     plot(Y{k}(1,:),Y{k}(2,:),'k'); hold on;
    plot(X_cond{k}{2}(1,2:end),X_cond{k}{2}(2,2:end),'r'); hold on;
end
title('Smoothed mode 2')
subplot(2,2,3);
for k=Tracks
     plot(Y{k}(1,:),Y{k}(2,:),'k'); hold on;
    plot(X_cond{k}{3}(1,2:end),X_cond{k}{3}(2,2:end),'r'); hold on;
end
title('Smoothed mode 3')
subplot(2,2,4);
for k=Tracks
     plot(Y{k}(1,:),Y{k}(2,:),'k'); hold on;
    plot(X_out{k}(1,2:end),X_out{k}(2,2:end),'r'); hold on;
end
title('Smoothed merged')

%% Convert mask into a surface 
AA = double(BW); % create a surface
Xx = 1:1:size(A,2); % create the grid of x coords
Yy = 1:1:size(A,1); % create the grid of y coords
[Yy_grid,Xx_grid] = meshgrid(Xx,Yy); % mesh
white=[1,1,1]; % surface colour
gg = [0.8,0.8,0.8]; % extra colour for cells
%%
% fig3
% Cell tracks colored with modes
fig('Modes',visFlag)
imshow(A); hold on;
% hold on;
for j = Tracks
   plot(X_out{j}(1,:),X_out{j}(2,:),'-b','LineWidth',1); hold on; 
   tx = num2str(j);
   text(X_out{j}(1,1),X_out{j}(2,1),2, tx, 'interpreter', 'latex','Color','k','FontSize',12);
end
 for j = 1:nResampled
   plot(X_resample{j}(1:end-1,1),X_resample{j}(1:end-1,2),'-g','LineWidth',1); hold on; 
 end
  for j = 1:nDead
   plot(X_dead{j}(1:end-1,1),X_dead{j}(1:end-1,2),'-','Color',white,'LineWidth',1); hold on; 
  end
xlim(x_lim);ylim(y_lim);
hold on;
surf(Yy_grid,Xx_grid,-AA,'FaceColor',white,'EdgeColor',white);
view(2)
hold on;
hold on;
text(625+ padW, 825 + padH,2, '$\star$  ', 'interpreter', 'latex','Color','g','FontSize',20);
text(650+ padW, 825 + padH,2, ' - $M^1$', 'interpreter', 'latex','Color','k','FontSize',20);
text(750+ padW, 825 + padH,2, '$\star$  ', 'interpreter', 'latex','Color','b','FontSize',20);
text(775+ padW, 825 + padH,2, ' - $M^2$', 'interpreter', 'latex','Color','k','FontSize',20);
text(875+ padW, 825 + padH,2, '$\star$  ', 'interpreter', 'latex','Color',white,'FontSize',20);
text(900+ padW, 825 + padH,2, ' - $M^3$', 'interpreter', 'latex','Color','k','FontSize',20);
set(gca,'Ydir','reverse')
% line([50,221],[825,825],'Color','k','LineWidth',5);
 line([250,250+100*cc],[y_max-20,y_max-20],[2,2],'Color','k','LineWidth',5);
% line([50,221],[575,575],'Color','k','LineWidth',5);
txt = ('100 $\mu$m');
text(250,y_max-45, 2,txt,'Color','k','FontSize',20)
text(70,550,2,txt,'Color','k','FontSize',20)
set(gca,'Ydir','reverse')
print([FigFolder,'modes_all_',Injury,num2str(iFish)],saveFormat)
%% Cell tracks colored with modes
% Mode 1
fig('Mode 1',visFlag);
imshow(A); hold on;
for k = Tracks
   c = Mu_out{k}(1,:);
    cplot(X_out{k}(1,:),X_out{k}(2,:),c,'-',...
    'markerfacecolor','flat',...
    'markeredgecolor','w','linewidth',1.0); hold on;
end
colorbar;
caxis([0 1]);
xlim(x_lim);ylim(y_lim);
hold on;
surf(Yy_grid,Xx_grid,-AA,'FaceColor',white,'EdgeColor',white);
view(2)
hold on;
line([250,250+100*cc],[y_max-20,y_max-20],[2,2],'Color','k','LineWidth',5);
txt = ('100 $\mu$m');
text(250,y_max-45, 2,txt,'Color','k','FontSize',20)
set(gca,'Ydir','reverse')
tightfig;
print([FigFolder,'mode_1_',Injury,num2str(iFish)],saveFormat)

% Mode 2
fig('Mode 2',visFlag);
imshow(A); hold on;
for k = Tracks
   c = Mu_out{k}(2,:);
    cplot(X_out{k}(1,:),X_out{k}(2,:),c,'-',...
    'markerfacecolor','flat',...
    'markeredgecolor','w','linewidth',1.0); hold on;
end
colorbar;
caxis([0 1]);
xlim(x_lim);ylim(y_lim);
hold on;
surf(Yy_grid,Xx_grid,-AA,'FaceColor',white,'EdgeColor',white);
view(2)
hold on;
line([250,250+100*cc],[y_max-20,y_max-20],[2,2],'Color','k','LineWidth',5);
txt = ('100 $\mu$m');
text(250,y_max-45, 2,txt,'Color','k','FontSize',20)
set(gca,'Ydir','reverse')
tightfig;
print([FigFolder,'mode_2_',Injury,num2str(iFish)],saveFormat)

% Mode 3
fig('Mode 3',visFlag);
imshow(A); hold on;
for k = Tracks
   c = Mu_out{k}(3,:);
    cplot(X_out{k}(1,:),X_out{k}(2,:),c,'-',...
    'markerfacecolor','flat',...
    'markeredgecolor','w','linewidth',1.0); hold on;
end
colorbar;
caxis([0 1]);
xlim(x_lim);ylim(y_lim);
hold on;
surf(Yy_grid,Xx_grid,-AA,'FaceColor',white,'EdgeColor',white);
view(2)
hold on;
 line([250,250+100*cc],[y_max-20,y_max-20],[2,2],'Color','k','LineWidth',5);
txt = ('100 $\mu$m');
text(250,y_max-45, 2,txt,'Color','k','FontSize',24)
set(gca,'Ydir','reverse')
tightfig;
print([FigFolder,'mode_3_',Injury,num2str(iFish)],saveFormat)
%% fig6 - Heatmap
fig('Heatmap',visFlag);
imshow(A); hold on;
plot_heatmap(Theta,Z,knots,grid_limits,basis_type);
% alpha(0.5)
hold on;
surf(Yy_grid,Xx_grid,-AA,'FaceColor',white,'EdgeColor',white);
view(2)
colorbar;
hold on;
xlim(x_lim);ylim(y_lim);
line([250,250+100*cc],[y_max-20,y_max-20],[2,2],'Color','k','LineWidth',5);
% line([50,221],[575,575],'Color','k','LineWidth',5);
txt = ('100 $\mu$m');
text(250,y_max-45, 2,txt,'Color','k','FontSize',20)
% text(70,550, 2,txt,'Color','k','FontSize',20)
set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(gca,'Ydir','reverse');
print([FigFolder,'heatmap_',Injury,num2str(iFish)],saveFormat)
