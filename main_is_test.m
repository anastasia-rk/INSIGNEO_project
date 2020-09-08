local_init;
%% Select the pattern
nTracks = 100;                                                              % number of tracks to generate in each map
% pattern = questdlg('Starting position distribution', ...
%     'Select pattern',...
% 	'Uniform','Line','Point','');
pattern = 'Uniform';
switch pattern
    case 'Uniform'
         folderName = 'Simulated/uniform_start/';
    case 'Line'
         folderName = 'Simulated/line_start/';
    case 'Point'
         folderName = 'Simulated/point_start/';
end
%% Setup parameters
%basis_type = 'gaussian';
basis_type = 'bspline';
% Set up limits of the grid: x_min,y_min,x_max,y_max
grid_limits = [0, 0, 1000, 1000];
% Set up number of basis functions
nx = 4; ny = 4; order = 4;
[knots] = setup_spline_support(grid_limits,nx,ny,order);
Z = 0;
ll = size(knots,2)/2;
%% Initialise field estimate and model parameters
grid_limits1 = [100, 200, 720, 750];
Theta = initiate_field(0,100,grid_limits1,knots);                            %fit splines to a flat surface
Theta_old = Theta;

 Theta_model = ones(ll,1);
 Theta = ones(ll,1);
 Theta_model(1:4,1) = 10;
 Theta_model(5:8,1) = 100;
 Theta_model(9:12,1) = 190;
 Theta_model(13:16,1) = 290;
 mu_field = 1;
 %% State-space model parameters
x_len   = 4;
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
C = [I O];
% Disturbance matrix
G_cv = [T^2/2*I; T*I];
G_rw = [T^2/2*I; T*I]; 
%% Covariances for the Kalman filters
sig1_Q = 0.5; % RW  disturbance - cell speed
sig2_Q = 0.2; % CV disturbance - random component of cell acceleration
% Process noise - for random walk
mm = sig1_Q*G_rw*G_rw';  
dd = diag(mm);
Q_rw = diag(dd);
Q_rw_inv = (1/sig1_Q)*eye(2);
% Process noise - for responsive mode
mm = sig2_Q*G_cv*G_cv';
dd = diag(mm);
Q_cv = diag(dd);
Q_cv_inv = (1/sig2_Q)*eye(2);
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
mu_0 = ones(1,n_models);
mu_0 = mu_0./n_models;
% for the maximisation step
Sig_w{1} = G{1}*(inv(G{1}'*G{1}))'*Q_cv_inv*inv(G{1}'*G{1})*G{1}';
Sig_w{2} = G{2}*(inv(G{2}'*G{2}))'*Q_rw_inv*inv(G{2}'*G{2})*G{2}';
%% Monte Carlo loop
L = 20;                                                                     % number of particles
models = [1];                                                               % Monte-Carlo simulations
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local');
end
for iModel = models
clear X Y Mode_model
load([folderName,'simulated_tracks_' num2str(iModel)]);
Tracks = [1:nTracks]; 
for k=Tracks % Assign track lengths outside of parfor loop for use in sliced variables
    T(k) = length(Y{k});
end
converged_theta = false;
converged_phi   = false;
converged_r     = false;
iter_max  = 10;
iter      = 0;
Theta = zeros(ll,1); 
p_tr = 0.5*ones(n_models);
%% EM loop
while  (iter < iter_max)
%% Expectation step
fprintf('E-step... \n')
iter = iter + 1;  
%% Parfor loop for IMM
tic
for k=Tracks     
%% recurtion cycle for IMM filter
    [x_m,P,mu,mode_probable_f{k}] = IMM_forward(T(k),x_len,n_models,p_tr,mu_0,Y{k},Q,R,F,B,C,Theta,knots,order);
%% Recurtion cycle for IMM RTS smoother
    [x_s,P_s,mu_s,x_merged,P_merged,mu_cond{k},mode_probable{k}] = IMM_backward(T(k),x_len,n_models,p_tr,mu,x_m,P,F,B,Q,Theta,knots,order);
%% Mixture elements
    [mixture{k}] = gm_elements(x_s,P_s,mu_s);
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
        phi_est{j,l}(iter,iModel) = transitions{iter}(j,l);
    end % for modes (l)
end % for modes (j)
p_tr_old = p_tr;
p_tr = transitions{iter}
%% Field parameters
sum1 = 0; sum2 =  0; sum3 = 0;
tic
parfor k=Tracks
     for t=1:T(k)-1
         % sample the particles from proposal pdf
          x_tilde = mvnrnd(X_out{k}(:,t)',P_out{k}{t},L);
%          [x_tilde,denoms] = univmultivar(X_out{k}(:,t),P_out{k}{t},L);
         % compute values of proposal for each particle
%          denoms = mvnpdf(x_tilde,X_out{k}(:,t)',P_out{k}{t});
         % compute true pdf for each particle
%           numers = pdf(mixture{k}{t},x_tilde);
%          numers = mvnpdf(x_tilde,X_out{k}(:,t)',P_out{k}{t});
         % compute particle weights
         weights =  ones(size(numers));  % 
%          weights = weights/sum(weights);
         % transpose the vector to use in the dynamical function readily
         x_tilde = x_tilde'; 
         grads = gradient_compute(x_tilde,knots,order);
         % evaluate the gradient at each particle
         for j=1:n_models
            % sample the particles from proposal pdf
             x_j_tilde = mvnrnd(X_cond{k}{j}(:,t+1)',P_cond{k}{j,t+1},L);
%             [x_j_tilde,denoms_j] = univmultivar(X_cond{k}{j}(:,t+1),P_cond{k}{j,t+1},L)
            % compute values of proposal for each particle
%             denoms_j = mvnpdf(x_j_tilde,X_cond{k}{j}(:,t+1)',P_cond{k}{j,t+1});
            % compute true pdf for each particle
%             numers_j = pdf(mixture{k}{t+1},x_j_tilde);
%             numers_j = mvnpdf(x_j_tilde,X_cond{k}{j}(:,t+1)',P_cond{k}{j,t+1});
            % compute particle weights
            weights_j =  ones(size(numers_j)); % numers_j./1; % 
%             weights_j = weights_j/sum(weights_j);
            % transpose the vector to use in the dynamical function readily
            x_j_tilde = x_j_tilde'; 
            for i=1:L
                bb = B{j}*mu_field* grads{i};
                sum1 = sum1 + weights(i)*Mu_out{k}(j,t)*bb'*Sig_w{j}*bb;
                sum3 = sum3 + weights(i)*Mu_out{k}(j,t)*bb'*Sig_w{j}*F{j}*x_tilde(:,i);
                for l=1:L
                    sum2 = sum2 + weights(i)*weights_j(l)*Mu_out{k}(j,t)*bb'*Sig_w{j}*x_j_tilde(:,l);
                end % for particles at x_tilde_j (l)
            end  % for particles at x_tilse (i)
         end % for modes (j) 
     end % for times (t)
end % for tracks (k)
toc
%% Normalise integrals using weights
sum1 = sum1/L;
sum2 = sum2/(L^2);
sum3 = sum3/L;    
% Parameter MLE
Theta_old = Theta;
Theta = pinv(sum1)*(sum2-sum3);
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
%% Save MC iteration results
fprintf([num2str(iModel) ' done... \n'])
% Save model paremeters
Theta_all(:,iModel) = Theta;
Mode_all{iModel} = mode_probable;
Model_model_all{iModel} = Mode_model;
Fisher_info{iModel} = Expected_info;
Miss_fisher{iModel} = Missing_info;
end % for MC simulation (iModel)
delete(pool);
%%
fiName = ['Results/is_',pattern];
save(fiName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example track - filtered, smoothed, conditioned
fig('Track',visFlag); 
subplot(3,1,1);
plot(x_plot(1,2:end),x_plot(2,2:end)); hold on;
plot(x_sm_plot(1,:),x_sm_plot(2,:)); hold on;
plot(X_cond{100,1}(1,:),X_cond{100,1}(2,:)); hold on;
plot(X_cond{100,2}(1,:),X_cond{100,2}(2,:)); hold on;
subplot(3,1,2);
plot(x_plot(3,:)); hold on;
plot(x_sm_plot(3,:)); hold on;
plot(X_cond{100,1}(3,:)); hold on;
plot(X_cond{100,2}(3,:)); hold on;
subplot(3,1,3); 
plot(x_plot(4,:)); hold on;
plot(x_sm_plot(4,:)); hold on;
plot(X_cond{100,1}(4,:)); hold on;
plot(X_cond{100,2}(4,:)); hold on;
%% Transition probability matrix for all models - statistics
iter_plot = iter+1;
for j=1:n_models
    for l=1:n_models
        phi_mean{j,l}(1) = 0.5;
        for it = 1:iter_plot
            vec = nonzeros(phi_est{j,l}(it,:));
            phi_mean{j,l}(it+1) = mean(vec);
            phi_std{j,l}(it+1) = std(vec);  
        end
    end
end
fig('Convergence of Phi',visFlag); 
figure; 
% line([0, iter_plot],[0.9, 0.9],'Color','black'); hold on;
% line([0, iter_plot],[0.8, 0.8],'Color','black'); hold on;
% line([0, iter_plot],[0.2, 0.2],'Color','black'); hold on;
% line([0, iter_plot],[0.1, 0.1],'Color','black'); hold on;
x_patch = [0:iter_plot];
xxx = [x_patch fliplr(x_patch)];
iPhi = 0;
for j=1:n_models
    for l=1:n_models
        iPhi = iPhi + 1;
        thisline(iPhi) = plot(x_patch,phi_mean{j,l}); hold on;
        color = get(thisline(iPhi), 'Color');
        yy = [phi_mean{j,l}+3*phi_std{j,l} fliplr(phi_mean{j,l}-3*phi_std{j,l})];
        patch(xxx,yy,color,'EdgeColor',color,'FaceAlpha',0.2); hold on;
%         for it=1:iter_plot
%             scatter(it*ones(1,iModel),phi_est{j,l}(it,:),'filled','MarkerFaceColor',color); hold on %
%         end
        clear color
    end
end
ylim([0,1]);
xlabel('Iteration of the EM algorithm');
ylabel('$\hat{\phi}_{j,i}$')
legend(thisline,{'$\hat{\phi}_{1,1}$','$\hat{\phi}_{1,2}$','$\hat{\phi}_{2,1}$','$\hat{\phi}_{2,2}$'})
print([FigFolder,'phi_convergence_',pattern],saveFormat)
tikzName = [TikzFolder,'transitions.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '10cm', 'width','10cm','checkForUpdates',false);
%% Modes
k = 38
fig('Mode for 1',visFlag); 
subplot(2,1,1)
plot(Mode_model{k},'','LineWidth',2); hold on;
plot(mode_probable{k},'--','Linewidth',2);
ylim([0 3]);
legend('Modelled mode','Estimated mode');
ylabel('$\mathbf{r}^{k}_{t}$'); xlabel('Time, min')
subplot(2,1,2)
plot(Mu_out{k}(1,:)); hold on;
plot(Mu_out{k}(2,:));
legend('j =  1','j = 2')
ylabel('$\mu^{k}_{t}(j)$'); xlabel('Time, min')
tikzName = [TikzFolder,'single_agent_mode.tikz']; 
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '10cm', 'width','10cm','checkForUpdates',false);
 %% Mean estimate
for i=1:length(Theta)
    Theta_mean(i,1) = mean(Theta_all(i,:));
    Std_dev(i,1) = std(Theta_all(i,:));
end
for iModel=models
    Theta_temp = Theta_all(:,iModel);
    Dist(iModel) = norm(Theta_model - Theta_temp);
end
[D_min,I_min] = min(Dist);
Theta_closest =  Theta_all(:,I_min);
[D_min,I_max] = max(Dist);
Theta_furthest =  Theta_all(:,I_max);
Theta_median = median(Theta_all');
%% Plot absolute bias
Theta_diff = Theta_model - Theta_mean;
fig('Bias',visFlag); 
colormap(my_map);
plot_heatmap(Theta_diff,Z,knots,grid_limits,basis_type);
hold on;
side = knots(1,2) - knots(1,1);
for i=1:length(Theta_mean)
    a = knots(1,i*2-1);
    b = knots(2,i*2-1);
    p=rectangle('Position',[a,b,side,side],'Curvature',0.1,'EdgeColor','w'); 
    hold on;
end
xlabel('\textrm{X, a.u.}', 'interpreter', 'latex');
ylabel('\textrm{Y,  a.u.}', 'interpreter', 'latex');
print([FigFolder,'bias_',pattern],saveFormat)
%% Plot estimated field
fig('Field est.',visFlag); 
colormap(my_map);
plot_heatmap(Theta_mean,Z,knots,grid_limits,basis_type);
hold on;
side = knots(1,2) - knots(1,1);
for i=1:length(Theta)
    a = knots(1,i*2-1);
    b = knots(2,i*2-1);
    p=rectangle('Position',[a,b,side,side],'Curvature',0.1,'EdgeColor','w'); 
    hold on;
end
xlabel('\textrm{X, a.u.}', 'interpreter', 'latex');
ylabel('\textrm{Y,  a.u.}', 'interpreter', 'latex');
print([FigFolder,'mean_estimate_',pattern],saveFormat)
%% Plot model gradient
fig('Grad mod.',visFlag); 
colormap(my_map);
plot_gradient(Theta_model,Z,knots,grid_limits,basis_type);
hold on;
switch pattern
    case 'Uniform'
        hold off;
    case 'Line'
         line([450, 450],[450, 450],'Color','m','LineWidth',4);
    case 'Point'
         line([450, 450],[250, 750],'Color','m','LineWidth',4);
end
colorbar
xlim([0 1000]);
ylim([0 1000]);
xlabel('\textrm{X,  a.u.}', 'interpreter', 'latex');
ylabel('\textrm{Y,  a.u.}', 'interpreter', 'latex');
print([FigFolder,'model_gradient'],saveFormat)
%% Plot estimated gradient
fig('Grad est.',visFlag); 
colormap(my_map);
plot_gradient(Theta_mean,Z,knots,grid_limits,basis_type);local_init;
%% Select the pattern
nTracks = 100;                                                              % number of tracks to generate in each map
% pattern = questdlg('Starting position distribution', ...
%     'Select pattern',...
% 	'Uniform','Line','Point','');
pattern = 'Uniform';
switch pattern
    case 'Uniform'
         folderName = 'Simulated/uniform_start/';
    case 'Line'
         folderName = 'Simulated/line_start/';
    case 'Point'
         folderName = 'Simulated/point_start/';
end
%% Setup parameters
%basis_type = 'gaussian';
basis_type = 'bspline';
% Set up limits of the grid: x_min,y_min,x_max,y_max
grid_limits = [0, 0, 1000, 1000];
% Set up number of basis functions
nx = 4; ny = 4; order = 4;
[knots] = setup_spline_support(grid_limits,nx,ny,order);
Z = 0;
ll = size(knots,2)/2;
%% Initialise field estimate and model parameters
grid_limits1 = [100, 200, 720, 750];
Theta = initiate_field(0,100,grid_limits1,knots);                            %fit splines to a flat surface
Theta_old = Theta;

 Theta_model = ones(ll,1);
 Theta = ones(ll,1);
 Theta_model(1:4,1) = 10;
 Theta_model(5:8,1) = 100;
 Theta_model(9:12,1) = 190;
 Theta_model(13:16,1) = 290;
 mu_field = 1;
 %% State-space model parameters
x_len   = 4;
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
C = [I O];
% Disturbance matrix
G_cv = [T^2/2*I; T*I];
G_rw = [T^2/2*I; T*I]; 
%% Covariances for the Kalman filters
sig1_Q = 0.5; % RW  disturbance - cell speed
sig2_Q = 0.2; % CV disturbance - random component of cell acceleration
% Process noise - for random walk
mm = sig1_Q*G_rw*G_rw';  
dd = diag(mm);
Q_rw = diag(dd);
Q_rw_inv = (1/sig1_Q)*eye(2);
% Process noise - for responsive mode
mm = sig2_Q*G_cv*G_cv';
dd = diag(mm);
Q_cv = diag(dd);
Q_cv_inv = (1/sig2_Q)*eye(2);
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
mu_0 = ones(1,n_models);
mu_0 = mu_0./n_models;
% for the maximisation step
Sig_w{1} = G{1}*(inv(G{1}'*G{1}))'*Q_cv_inv*inv(G{1}'*G{1})*G{1}';
Sig_w{2} = G{2}*(inv(G{2}'*G{2}))'*Q_rw_inv*inv(G{2}'*G{2})*G{2}';
%% Monte Carlo loop
L = 20;                                                                     % number of particles
models = [1];                                                            % Monte-Carlo simulations
pool =  parpool('local');
for iModel = models
clear X Y Mode_model
load([folderName,'simulated_tracks_' num2str(iModel)]);
Tracks = [1:nTracks]; 
for k=Tracks % Assign track lengths outside of parfor loop for use in sliced variables
    T(k) = length(Y{k});
end
converged_theta = false;
converged_phi   = false;
converged_r     = false;
iter_max  = 10;
iter      = 0;
Theta = zeros(ll,1); 
p_tr = 0.5*ones(n_models);
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
%% Mixture elements
    [mixture{k}] = gm_elements(x_s,P_s,mu_s);
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
        phi_est{j,l}(iter,iModel) = transitions{iter}(j,l);
    end % for modes (l)
end % for modes (j)
p_tr_old = p_tr;
p_tr = transitions{iter}
%% Field parameters
sum1 = 0; sum2 =  0; sum3 = 0;
tic
parfor k=Tracks
     for t=1:T(k)-1
         % sample the particles from proposal pdf
         x_tilde = mvnrnd(X_out{k}(:,t)',P_out{k}{t},L);
         % compute values of proposal for each particle
         denoms = mvnpdf(x_tilde,X_out{k}(:,t)',P_out{k}{t});
         % compute true pdf for each particle
         numers = pdf(mixture{k}{t},x_tilde);
         % compute particle weights
         weights = numers./denoms;
         % transpose the vector to use in the dynamical function readily
         x_tilde = x_tilde'; 
         grads = gradient_compute(x_tilde,knots,order);
         % evaluate the gradient at each particle
         for j=1:n_models
            % sample the particles from proposal pdf
            x_j_tilde = mvnrnd(X_cond{k}{j}(:,t+1)',P_cond{k}{j,t+1},L);
            % compute values of proposal for each particle
            denoms_j = mvnpdf(x_j_tilde,X_cond{k}{j}(:,t+1)',P_cond{k}{j,t+1});
            % compute true pdf for each particle
            numers_j = pdf(mixture{k}{t+1},x_j_tilde);
            % compute particle weights
            weights_j = numers_j./denoms_j;
            % transpose the vector to use in the dynamical function readily
            x_j_tilde = x_j_tilde'; 
            for i=1:L
                bb = B{j}*mu_field* grads{i};
                sum1 = sum1 + weights(i)*Mu_out{k}(j,t)*bb'*Sig_w{j}*bb;
                sum3 = sum3 + weights(i)*Mu_out{k}(j,t)*bb'*Sig_w{j}*F{j}*x_tilde(:,i);
		for l=1:L
                    sum2 = sum2 + weights(i)*weights_j(l)*Mu_out{k}(j,t)*bb'*Sig_w{j}*x_j_tilde(:,l);
                end % for particles at x_tilde_j (l)
            end  % for particles at x_tilse (i)
         end % for modes (j) 
     end % for times (t)
end % for tracks (k)
toc
%% Normalise integrals using weights
sum1 = sum1/L;
sum2 = sum2/(L^2);
sum3 = sum3/L;    
% Parameter MLE
Theta_old = Theta;
Theta = pinv(sum1)*(sum2-sum3)
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
%% Save MC iteration results
fprintf([num2str(iModel) ' done... \n'])
% Save model paremeters
Theta_all(:,iModel) = Theta;
Mode_all{iModel} = mode_probable;
Model_model_all{iModel} = Mode_model;
Fisher_info{iModel} = Expected_info;
Miss_fisher{iModel} = Missing_info;
end % for MC simulation (iModel)
delete(pool);
fiName = ['Results/is_',pattern];
save(finName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example track - filtered, smoothed, conditioned
fig('Track',visFlag); 
subplot(3,1,1);
plot(x_plot(1,2:end),x_plot(2,2:end)); hold on;
plot(x_sm_plot(1,:),x_sm_plot(2,:)); hold on;
plot(X_cond{100,1}(1,:),X_cond{100,1}(2,:)); hold on;
plot(X_cond{100,2}(1,:),X_cond{100,2}(2,:)); hold on;
subplot(3,1,2);
plot(x_plot(3,:)); hold on;
plot(x_sm_plot(3,:)); hold on;
plot(X_cond{100,1}(3,:)); hold on;
plot(X_cond{100,2}(3,:)); hold on;
subplot(3,1,3); 
plot(x_plot(4,:)); hold on;
plot(x_sm_plot(4,:)); hold on;
plot(X_cond{100,1}(4,:)); hold on;
plot(X_cond{100,2}(4,:)); hold on;
%% Transition probability matrix for all models - statistics
iter_plot = iter+1;
for j=1:n_models
    for l=1:n_models
        phi_mean{j,l}(1) = 0.5;
        for it = 1:iter_plot
            vec = nonzeros(phi_est{j,l}(it,:));
            phi_mean{j,l}(it+1) = mean(vec);
            phi_std{j,l}(it+1) = std(vec);  
        end
    end
end
fig('Convergence of Phi',visFlag); 
figure; 
% line([0, iter_plot],[0.9, 0.9],'Color','black'); hold on;
% line([0, iter_plot],[0.8, 0.8],'Color','black'); hold on;
% line([0, iter_plot],[0.2, 0.2],'Color','black'); hold on;
% line([0, iter_plot],[0.1, 0.1],'Color','black'); hold on;
x_patch = [0:iter_plot];
xxx = [x_patch fliplr(x_patch)];
iPhi = 0;
for j=1:n_models
    for l=1:n_models
        iPhi = iPhi + 1;
        thisline(iPhi) = plot(x_patch,phi_mean{j,l}); hold on;
        color = get(thisline(iPhi), 'Color');
        yy = [phi_mean{j,l}+3*phi_std{j,l} fliplr(phi_mean{j,l}-3*phi_std{j,l})];
        patch(xxx,yy,color,'EdgeColor',color,'FaceAlpha',0.2); hold on;
%         for it=1:iter_plot
%             scatter(it*ones(1,iModel),phi_est{j,l}(it,:),'filled','MarkerFaceColor',color); hold on %
%         end
        clear color
    end
end
ylim([0,1]);
xlabel('Iteration of the EM algorithm');
ylabel('$\hat{\phi}_{j,i}$')
legend(thisline,{'$\hat{\phi}_{1,1}$','$\hat{\phi}_{1,2}$','$\hat{\phi}_{2,1}$','$\hat{\phi}_{2,2}$'})
print([FigFolder,'phi_convergence_',pattern],saveFormat)
tikzName = [TikzFolder,'transitions.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '10cm', 'width','10cm','checkForUpdates',false);
%% Modes
k = 38
fig('Mode for 1',visFlag); 
subplot(2,1,1)
plot(Mode_model{k},'','LineWidth',2); hold on;
plot(mode_probable{k},'--','Linewidth',2);
ylim([0 3]);
legend('Modelled mode','Estimated mode');
ylabel('$\mathbf{r}^{k}_{t}$'); xlabel('Time, min')
subplot(2,1,2)
plot(Mu_out{k}(1,:)); hold on;
plot(Mu_out{k}(2,:));
legend('j =  1','j = 2')
ylabel('$\mu^{k}_{t}(j)$'); xlabel('Time, min')
tikzName = [TikzFolder,'single_agent_mode.tikz']; 
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '10cm', 'width','10cm','checkForUpdates',false);
 %% Mean estimate
for i=1:length(Theta)
    Theta_mean(i,1) = mean(Theta_all(i,:));
    Std_dev(i,1) = std(Theta_all(i,:));
end
for iModel=models
    Theta_temp = Theta_all(:,iModel);
    Dist(iModel) = norm(Theta_model - Theta_temp);
end
[D_min,I_min] = min(Dist);
Theta_closest =  Theta_all(:,I_min);
[D_min,I_max] = max(Dist);
Theta_furthest =  Theta_all(:,I_max);
Theta_median = median(Theta_all');
%% Plot absolute bias
Theta_diff = Theta_model - Theta_mean;
fig('Bias',visFlag); 
colormap(my_map);
plot_heatmap(Theta_diff,Z,knots,grid_limits,basis_type);
hold on;
side = knots(1,2) - knots(1,1);
for i=1:length(Theta_mean)
    a = knots(1,i*2-1);
    b = knots(2,i*2-1);
    p=rectangle('Position',[a,b,side,side],'Curvature',0.1,'EdgeColor','w'); 
    hold on;
end
xlabel('\textrm{X, a.u.}', 'interpreter', 'latex');
ylabel('\textrm{Y,  a.u.}', 'interpreter', 'latex');
print([FigFolder,'bias_',pattern],saveFormat)
%% Plot estimated field
fig('Field est.',visFlag); 
colormap(my_map);
plot_heatmap(Theta_mean,Z,knots,grid_limits,basis_type);
hold on;
side = knots(1,2) - knots(1,1);
for i=1:length(Theta)
    a = knots(1,i*2-1);
    b = knots(2,i*2-1);
    p=rectangle('Position',[a,b,side,side],'Curvature',0.1,'EdgeColor','w'); 
    hold on;
end
xlabel('\textrm{X, a.u.}', 'interpreter', 'latex');
ylabel('\textrm{Y,  a.u.}', 'interpreter', 'latex');
print([FigFolder,'mean_estimate_',pattern],saveFormat)
%% Plot model gradient
fig('Grad mod.',visFlag); 
colormap(my_map);
plot_gradient(Theta_model,Z,knots,grid_limits,basis_type);
hold on;
switch pattern
    case 'Uniform'
        hold off;
    case 'Line'
         line([450, 450],[450, 450],'Color','m','LineWidth',4);
    case 'Point'
         line([450, 450],[250, 750],'Color','m','LineWidth',4);
end
colorbar
xlim([0 1000]);
ylim([0 1000]);
xlabel('\textrm{X,  a.u.}', 'interpreter', 'latex');
ylabel('\textrm{Y,  a.u.}', 'interpreter', 'latex');
print([FigFolder,'model_gradient'],saveFormat)
%% Plot estimated gradient
fig('Grad est.',visFlag); 
colormap(my_map);
plot_gradient(Theta_mean,Z,knots,grid_limits,basis_type);
hold on;
switch pattern
    case 'Uniform'
        hold off;
    case 'Line'
         line([450, 450],[450, 450],'Color','m','LineWidth',4);
    case 'Point'
         line([450, 450],[250, 750],'Color','m','LineWidth',4);
end
colorbar
xlim([0 1000]);
ylim([0 1000]);
xlabel('\textrm{X,  a.u.}', 'interpreter', 'latex');
ylabel('\textrm{Y,  a.u.}', 'interpreter', 'latex');
print([FigFolder,'estimated_gradient',pattern],saveFormat)
%% For the final MC realisation - modelled tracks and modes
fig('Modes mod.',visFlag); 
colormap(my_map);
plot_gradient(Theta,Z,knots,grid_limits,basis_type);
hold on;
for k=1:20
    c = Mode_model{k}(1,1:end-1);
    cplot(X{k}(1,:),X{k}(2,:),c,'-',...
    'markerfacecolor','flat',...
    'markeredgecolor','w','linewidth',2.0); hold on;
end
colorbar;
caxis([0 2]);
xlim([230 700]);
ylim([300 800]);
xlabel('\textrm{X,  a.u.}');
ylabel('\textrm{Y,  a.u.}');
print([FigFolder,'modes_modelled_',pattern],saveFormat)
tikzName = [TikzFolder,'modes_modelled',pattern,'.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '10cm', 'width','10cm','checkForUpdates',false);
%% For the final MC realisation - estimated tracks amd modes
fig('Modes est.',visFlag); 
colormap(my_map);
plot_gradient(Theta,Z,knots,grid_limits,basis_type);
hold on;
for k=1:20
    plot(X{k}(1,:),X{k}(2,:),'k'); hold on;
    c = Mu_out{k}(1,:);
    cplot(X_out{k}(1,:),X_out{k}(2,:),c,'-',...
    'markerfacecolor','flat',...
    'markeredgecolor','w','linewidth',2.0); hold on;
end
colorbar;
caxis([0 1]);
xlim([230 700]);
ylim([300 800]);
xlabel('\textrm{X,  a.u.}');
ylabel('\textrm{Y,  a.u.}');
print([FigFolder,'modes_estimated_',pattern],saveFormat)
tikzName = [TikzFolder,'modes_estimated_',pattern,'.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '10cm', 'width','10cm','checkForUpdates',false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variance ellipsoids (now 3D) - for thesis
% FIM100 = zeros(length(Theta),length(Theta));
% for iModel=models
%     FIM100 = FIM100 + Fisher_info{iModel};
% end
% COV100 = inv(FIM100);
% VOL100 = ellipsoid_volume(ll,0.05,COV100);
% TOTALV100 = trace(COV100)/length(Theta);
% %
% FIM50 = zeros(length(Theta),length(Theta));
% for iModel=1:50
%     FIM50 = FIM50 + Fisher_info{iModel};
% end
% COV50 = inv(FIM50);
% VOL50 = ellipsoid_volume(ll,0.05,COV50);
% for i=1:length(Theta)
%     Theta_mean_50(i,1) = mean(Theta_all(i,1:50));
% end
% TOTALV50 = trace(COV50)/length(Theta);
% %
% FIM25 = zeros(length(Theta),length(Theta));
% for iModel=1:25
%     FIM25 = FIM25 + Fisher_info{iModel};
% end
% COV25 = inv(FIM25);
% VOL25 = ellipsoid_volume(ll,0.05,COV25);
% for i=1:length(Theta)
%     Theta_mean_25(i,1) = mean(Theta_all(i,1:25));
% end
% TOTALV25 = trace(COV25)/length(Theta);
% %
% FIM10 = zeros(length(Theta),length(Theta));
% for iModel=1:10
%     FIM10 = FIM10 + Fisher_info{iModel};
% end
% COV10 = inv(FIM10);
% VOL10 = ellipsoid_volume(ll,0.05,COV10);
% for i=1:length(Theta)
%     Theta_mean_10(i,1) = mean(Theta_all(i,1:10));
% end
% TOTALV10 = trace(COV10)/length(Theta);
% %
% iModel = 1
% FIM1 = Fisher_info{iModel};
% COV1 = inv(FIM1);
% VOL1 = ellipsoid_volume(ll,0.05,COV1);
% for i=1:length(Theta)
%     Theta_mean_1(i,1) = mean(Theta_all(i,iModel));
% end
% TOTALV1 = trace(COV1)/length(Theta);
% %% Build ellipsoids for 100, 50, 25, 10 points
% i1 = 1;
% i2 = 2;
% i3 = 3;
%     Theta_hor = [Theta_all(i1,:);Theta_all(i2,:);Theta_all(i3,:)];
%     Theta_hor = Theta_hor';
%     covar(:,1) = COV100([i1,i2,i3],i1);
%     covar(:,2) = COV100([i1,i2,i3],i2);
%     covar(:,3) = COV100([i1,i2,i3],i3);
%     [X_el,Y_el,Z_el] = build_ellipsoid(Theta_hor,covar);
%     [X_el1,Y_el1,Z_el1] = build_ellipsoid_from_data(Theta_hor);
%     fig('100',visFlag); 
%     plot3(Theta_hor(:,1),Theta_hor(:,2),Theta_hor(:,3),'o','Color',[80/255, 80/255, 80/255]); hold on;
%     plot3(Theta_mean(i1),Theta_mean(i2),Theta_mean(i3),'*r','LineWidth',5); hold on;
%     plot3(Theta_model(i1),Theta_model(i2),Theta_model(i3),'*k','LineWidth',5); hold on;
%     h = surf(X_el,Y_el,Z_el); alpha 0.8   
%     set(h, 'facecolor',[243/255, 222/255, 187/255]);
%     set(h, 'edgecolor',[237/255, 177/255, 32/255]);
%     hold on;
%     h1 = surf(X_el1,Y_el1,Z_el1); alpha 0.5    
%     set(h1, 'facecolor',[187/255, 222/255, 243/255]);
%     set(h1, 'edgecolor',[32/255, 177/255, 237/255]);
%     grid on;
%     view(45, 25);
% %     xlim([-100,150]); ylim([-100,150]); zlim([-100,150]);
%     xlabel(['$\theta_{' num2str(i1) '}$']);
%     ylabel(['$\theta_{' num2str(i2) '}$']);
%     zlabel(['$\theta_{' num2str(i3) '}$']);
%     legend('Estimates','Mean estimate','True value','$95\%$ confidence region','$95\%$ covariance ellipsoid')
%     legend('Location','northeast')
% %%    
% Theta_hor = [Theta_all(i1,1:50);Theta_all(i2,1:50);Theta_all(i3,1:50)];
%     Theta_hor = Theta_hor';
%     covar(:,1) = COV50([i1,i2,i3],i1);
%     covar(:,2) = COV50([i1,i2,i3],i2);
%     covar(:,3) = COV50([i1,i2,i3],i3);
%     [X_el,Y_el,Z_el] = build_ellipsoid(Theta_hor,covar);
%     [X_el1,Y_el1,Z_el1] = build_ellipsoid_from_data(Theta_hor);?
%     fig('50',visFlag); 
%     plot3(Theta_hor(:,1),Theta_hor(:,2),Theta_hor(:,3),'o','Color',[80/255, 80/255, 80/255]); hold on;
%     plot3(Theta_mean_50(i1),Theta_mean_50(i2),Theta_mean_50(i3),'*r','LineWidth',5); hold on;
%     plot3(Theta_model(i1),Theta_model(i2),Theta_model(i3),'*k','LineWidth',5); hold on;
%     h = surf(X_el,Y_el,Z_el); alpha 0.6    
%     set(h, 'facecolor',[243/255, 222/255, 187/255]);
%     set(h, 'edgecolor',[237/255, 177/255, 32/255]);
%     hold on;
%     h1 = surf(X_el1,Y_el1,Z_el1); alpha 0.3    
%     set(h1, 'facecolor',[187/255, 222/255, 243/255]);
%     set(h1, 'edgecolor',[32/255, 177/255, 237/255]);
%     grid on;
%     view(45, 25);
% %     xlim([-100,150]); ylim([-100,150]); zlim([-100,150]);
%     xlabel(['$\theta_{' num2str(i1) '}$']);
%     ylabel(['$\theta_{' num2str(i2) '}$']);
%     zlabel(['$\theta_{' num2str(i3) '}$']);
%     legend('Estimates','Mean estimate','True value','$95\%$ confidence region','$95\%$ covariance ellipsoid')
%     legend('Location','northeast')% 
% %%    
% Theta_hor = [Theta_all(i1,1:25);Theta_all(i2,1:25);Theta_all(i3,1:25)];
%     Theta_hor = Theta_hor';
%     covar(:,1) = COV25([i1,i2,i3],i1);
%     covar(:,2) = COV25([i1,i2,i3],i2);
%     covar(:,3) = COV25([i1,i2,i3],i3);
%     [X_el,Y_el,Z_el] = build_ellipsoid(Theta_hor,covar);
%     [X_el1,Y_el1,Z_el1] = build_ellipsoid_from_data(Theta_hor);
%     fig('25',visFlag); 
%     plot3(Theta_hor(:,1),Theta_hor(:,2),Theta_hor(:,3),'o','Color',[80/255, 80/255, 80/255]); hold on;
%     plot3(Theta_mean_25(i1),Theta_mean_25(i2),Theta_mean_25(i3),'*r','LineWidth',5); hold on;
%     plot3(Theta_model(i1),Theta_model(i2),Theta_model(i3),'*k','LineWidth',5); hold on;
%     h = surf(X_el,Y_el,Z_el); alpha 0.8    
%     set(h, 'facecolor',[243/255, 222/255, 187/255]);
%     set(h, 'edgecolor',[237/255, 177/255, 32/255]);
%     hold on;
%     h1 = surf(X_el1,Y_el1,Z_el1); alpha 0.5    
%     set(h1, 'facecolor',[187/255, 222/255, 243/255]);
%     set(h1, 'edgecolor',[32/255, 177/255, 237/255]);
%     grid on;
%     view(45, 25);
% %     xlim([-100,150]); ylim([-100,150]); zlim([-100,150]);
%     xlabel(['$\theta_{' num2str(i1) '}$']);
%     ylabel(['$\theta_{' num2str(i2) '}$']);
%     zlabel(['$\theta_{' num2str(i3) '}$']);
%     legend('Estimates','Mean estimate','True value','$95\%$ confidence region','$95\%$ covariance ellipsoid')
%     legend('Location','northeast')
% %%    
% Theta_hor = [Theta_all(i1,1:10);Theta_all(i2,1:10);Theta_all(i3,1:10)];
%     Theta_hor = Theta_hor';
%     covar(:,1) = COV10([i1,i2,i3],i1);
%     covar(:,2) = COV10([i1,i2,i3],i2);
%     covar(:,3) = COV10([i1,i2,i3],i3);
%     [X_el,Y_el,Z_el] = build_ellipsoid(Theta_hor,covar);
%     [X_el1,Y_el1,Z_el1] = build_ellipsoid_from_data(Theta_hor);
%     fig('10',visFlag); 
%     plot3(Theta_hor(:,1),Theta_hor(:,2),Theta_hor(:,3),'o','Color',[80/255, 80/255, 80/255]); hold on;
%     plot3(Theta_mean_10(i1),Theta_mean_10(i2),Theta_mean_10(i3),'*r','LineWidth',5); hold on;
%     plot3(Theta_model(i1),Theta_model(i2),Theta_model(i3),'*k','LineWidth',5); hold on;
%     h = surf(X_el,Y_el,Z_el); alpha 0.8    
%     set(h, 'facecolor',[243/255, 222/255, 187/255]);
%     set(h, 'edgecolor',[237/255, 177/255, 32/255]);
%     hold on;
%     h1 = surf(X_el1,Y_el1,Z_el1); alpha 0.5    
%     set(h1, 'facecolor',[187/255, 222/255, 243/255]);
%     set(h1, 'edgecolor',[32/255, 177/255, 237/255]);
%     grid on;
%     view(45, 25);
% %     xlim([-100,150]); ylim([-100,150]); zlim([-100,150]);
%     xlabel(['$\theta_{' num2str(i1) '}$']);
%     ylabel(['$\theta_{' num2str(i2) '}$']);
%     zlabel(['$\theta_{' num2str(i3) '}$']);
%     legend('Estimates','Mean estimate','True value','$95\%$ confidence region','$95\%$ covariance ellipsoid')
%     legend('Location','northeast')
%% Pair-wise plots - uncomment if needed nore detail
% for iPair=1:2:ll-1
%     cc(:,1) = COV100([iPair, iPair+1],iPair);
%     cc(:,2) = COV100([iPair, iPair+1],iPair+1);
%     co{iPair} = cc;
%     mea{iPair} = [mean(Theta_all(iPair,:)),mean(Theta_all(iPair+1,:))];
%     ell_points_theta{iPair} = build_ellips(mea{iPair},co{iPair});
% end
% 
% for iPair=1:2:ll-1
%     fig(['Pair ',num2str(iPair)],visFlag); 
%     plot(Theta_all(iPair,:),Theta_all(iPair+1,:),'*k'); hold on;
%     plot(Theta_model(iPair),Theta_model(iPair+1),'*b','Linewidth',3); hold on;
%     plot(mea{iPair}(1),mea{iPair}(2),'*r','Linewidth',3); hold on;
%     plot(ell_points_theta{iPair}(:,1),ell_points_theta{iPair}(:,2),'r');
%     legend('MLEs','Mean MLE',)
%     xlabel(['$\Theta_{' num2str(iPair) '}$']);
%     ylabel(['$\Theta_{' num2str(iPair+1) '}$']);
% end

hold on;
switch pattern
    case 'Uniform'
        hold off;
    case 'Line'
         line([450, 450],[450, 450],'Color','m','LineWidth',4);
    case 'Point'
         line([450, 450],[250, 750],'Color','m','LineWidth',4);
end
colorbar
xlim([0 1000]);
ylim([0 1000]);
xlabel('\textrm{X,  a.u.}', 'interpreter', 'latex');
ylabel('\textrm{Y,  a.u.}', 'interpreter', 'latex');
print([FigFolder,'estimated_gradient',pattern],saveFormat)
%% For the final MC realisation - modelled tracks and modes
fig('Modes mod.',visFlag); 
colormap(my_map);
plot_gradient(Theta,Z,knots,grid_limits,basis_type);
hold on;
for k=1:20
    c = Mode_model{k}(1,1:end-1);
    cplot(X{k}(1,:),X{k}(2,:),c,'-',...
    'markerfacecolor','flat',...
    'markeredgecolor','w','linewidth',2.0); hold on;
end
colorbar;
caxis([0 2]);
xlim([230 700]);
ylim([300 800]);
xlabel('\textrm{X,  a.u.}');
ylabel('\textrm{Y,  a.u.}');
print([FigFolder,'modes_modelled_',pattern],saveFormat)
tikzName = [TikzFolder,'modes_modelled',pattern,'.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '10cm', 'width','10cm','checkForUpdates',false);
%% For the final MC realisation - estimated tracks amd modes
fig('Modes est.',visFlag); 
colormap(my_map);
plot_gradient(Theta,Z,knots,grid_limits,basis_type);
hold on;
for k=1:20
    plot(X{k}(1,:),X{k}(2,:),'k'); hold on;
    c = Mu_out{k}(1,:);
    cplot(X_out{k}(1,:),X_out{k}(2,:),c,'-',...
    'markerfacecolor','flat',...
    'markeredgecolor','w','linewidth',2.0); hold on;
end
colorbar;
caxis([0 1]);
xlim([230 700]);
ylim([300 800]);
xlabel('\textrm{X,  a.u.}');
ylabel('\textrm{Y,  a.u.}');
print([FigFolder,'modes_estimated_',pattern],saveFormat)
tikzName = [TikzFolder,'modes_estimated_',pattern,'.tikz'];
cleanfigure;
matlab2tikz(tikzName, 'showInfo', false,'parseStrings',false,'standalone', ...
            false, 'height', '10cm', 'width','10cm','checkForUpdates',false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variance ellipsoids (now 3D) - for thesis
% FIM100 = zeros(length(Theta),length(Theta));
% for iModel=models
%     FIM100 = FIM100 + Fisher_info{iModel};
% end
% COV100 = inv(FIM100);
% VOL100 = ellipsoid_volume(ll,0.05,COV100);
% TOTALV100 = trace(COV100)/length(Theta);
% %
% FIM50 = zeros(length(Theta),length(Theta));
% for iModel=1:50
%     FIM50 = FIM50 + Fisher_info{iModel};
% end
% COV50 = inv(FIM50);
% VOL50 = ellipsoid_volume(ll,0.05,COV50);
% for i=1:length(Theta)
%     Theta_mean_50(i,1) = mean(Theta_all(i,1:50));
% end
% TOTALV50 = trace(COV50)/length(Theta);
% %
% FIM25 = zeros(length(Theta),length(Theta));
% for iModel=1:25
%     FIM25 = FIM25 + Fisher_info{iModel};
% end
% COV25 = inv(FIM25);
% VOL25 = ellipsoid_volume(ll,0.05,COV25);
% for i=1:length(Theta)
%     Theta_mean_25(i,1) = mean(Theta_all(i,1:25));
% end
% TOTALV25 = trace(COV25)/length(Theta);
% %
% FIM10 = zeros(length(Theta),length(Theta));
% for iModel=1:10
%     FIM10 = FIM10 + Fisher_info{iModel};
% end
% COV10 = inv(FIM10);
% VOL10 = ellipsoid_volume(ll,0.05,COV10);
% for i=1:length(Theta)
%     Theta_mean_10(i,1) = mean(Theta_all(i,1:10));
% end
% TOTALV10 = trace(COV10)/length(Theta);
% %
% iModel = 1
% FIM1 = Fisher_info{iModel};
% COV1 = inv(FIM1);
% VOL1 = ellipsoid_volume(ll,0.05,COV1);
% for i=1:length(Theta)
%     Theta_mean_1(i,1) = mean(Theta_all(i,iModel));
% end
% TOTALV1 = trace(COV1)/length(Theta);
% %% Build ellipsoids for 100, 50, 25, 10 points
% i1 = 1;
% i2 = 2;
% i3 = 3;
%     Theta_hor = [Theta_all(i1,:);Theta_all(i2,:);Theta_all(i3,:)];
%     Theta_hor = Theta_hor';
%     covar(:,1) = COV100([i1,i2,i3],i1);
%     covar(:,2) = COV100([i1,i2,i3],i2);
%     covar(:,3) = COV100([i1,i2,i3],i3);
%     [X_el,Y_el,Z_el] = build_ellipsoid(Theta_hor,covar);
%     [X_el1,Y_el1,Z_el1] = build_ellipsoid_from_data(Theta_hor);
%     fig('100',visFlag); 
%     plot3(Theta_hor(:,1),Theta_hor(:,2),Theta_hor(:,3),'o','Color',[80/255, 80/255, 80/255]); hold on;
%     plot3(Theta_mean(i1),Theta_mean(i2),Theta_mean(i3),'*r','LineWidth',5); hold on;
%     plot3(Theta_model(i1),Theta_model(i2),Theta_model(i3),'*k','LineWidth',5); hold on;
%     h = surf(X_el,Y_el,Z_el); alpha 0.8   
%     set(h, 'facecolor',[243/255, 222/255, 187/255]);
%     set(h, 'edgecolor',[237/255, 177/255, 32/255]);
%     hold on;
%     h1 = surf(X_el1,Y_el1,Z_el1); alpha 0.5    
%     set(h1, 'facecolor',[187/255, 222/255, 243/255]);
%     set(h1, 'edgecolor',[32/255, 177/255, 237/255]);
%     grid on;
%     view(45, 25);
% %     xlim([-100,150]); ylim([-100,150]); zlim([-100,150]);
%     xlabel(['$\theta_{' num2str(i1) '}$']);
%     ylabel(['$\theta_{' num2str(i2) '}$']);
%     zlabel(['$\theta_{' num2str(i3) '}$']);
%     legend('Estimates','Mean estimate','True value','$95\%$ confidence region','$95\%$ covariance ellipsoid')
%     legend('Location','northeast')
% %%    
% Theta_hor = [Theta_all(i1,1:50);Theta_all(i2,1:50);Theta_all(i3,1:50)];
%     Theta_hor = Theta_hor';
%     covar(:,1) = COV50([i1,i2,i3],i1);
%     covar(:,2) = COV50([i1,i2,i3],i2);
%     covar(:,3) = COV50([i1,i2,i3],i3);
%     [X_el,Y_el,Z_el] = build_ellipsoid(Theta_hor,covar);
%     [X_el1,Y_el1,Z_el1] = build_ellipsoid_from_data(Theta_hor);?
%     fig('50',visFlag); 
%     plot3(Theta_hor(:,1),Theta_hor(:,2),Theta_hor(:,3),'o','Color',[80/255, 80/255, 80/255]); hold on;
%     plot3(Theta_mean_50(i1),Theta_mean_50(i2),Theta_mean_50(i3),'*r','LineWidth',5); hold on;
%     plot3(Theta_model(i1),Theta_model(i2),Theta_model(i3),'*k','LineWidth',5); hold on;
%     h = surf(X_el,Y_el,Z_el); alpha 0.6    
%     set(h, 'facecolor',[243/255, 222/255, 187/255]);
%     set(h, 'edgecolor',[237/255, 177/255, 32/255]);
%     hold on;
%     h1 = surf(X_el1,Y_el1,Z_el1); alpha 0.3    
%     set(h1, 'facecolor',[187/255, 222/255, 243/255]);
%     set(h1, 'edgecolor',[32/255, 177/255, 237/255]);
%     grid on;
%     view(45, 25);
% %     xlim([-100,150]); ylim([-100,150]); zlim([-100,150]);
%     xlabel(['$\theta_{' num2str(i1) '}$']);
%     ylabel(['$\theta_{' num2str(i2) '}$']);
%     zlabel(['$\theta_{' num2str(i3) '}$']);
%     legend('Estimates','Mean estimate','True value','$95\%$ confidence region','$95\%$ covariance ellipsoid')
%     legend('Location','northeast')% 
% %%    
% Theta_hor = [Theta_all(i1,1:25);Theta_all(i2,1:25);Theta_all(i3,1:25)];
%     Theta_hor = Theta_hor';
%     covar(:,1) = COV25([i1,i2,i3],i1);
%     covar(:,2) = COV25([i1,i2,i3],i2);
%     covar(:,3) = COV25([i1,i2,i3],i3);
%     [X_el,Y_el,Z_el] = build_ellipsoid(Theta_hor,covar);
%     [X_el1,Y_el1,Z_el1] = build_ellipsoid_from_data(Theta_hor);
%     fig('25',visFlag); 
%     plot3(Theta_hor(:,1),Theta_hor(:,2),Theta_hor(:,3),'o','Color',[80/255, 80/255, 80/255]); hold on;
%     plot3(Theta_mean_25(i1),Theta_mean_25(i2),Theta_mean_25(i3),'*r','LineWidth',5); hold on;
%     plot3(Theta_model(i1),Theta_model(i2),Theta_model(i3),'*k','LineWidth',5); hold on;
%     h = surf(X_el,Y_el,Z_el); alpha 0.8    
%     set(h, 'facecolor',[243/255, 222/255, 187/255]);
%     set(h, 'edgecolor',[237/255, 177/255, 32/255]);
%     hold on;
%     h1 = surf(X_el1,Y_el1,Z_el1); alpha 0.5    
%     set(h1, 'facecolor',[187/255, 222/255, 243/255]);
%     set(h1, 'edgecolor',[32/255, 177/255, 237/255]);
%     grid on;
%     view(45, 25);
% %     xlim([-100,150]); ylim([-100,150]); zlim([-100,150]);
%     xlabel(['$\theta_{' num2str(i1) '}$']);
%     ylabel(['$\theta_{' num2str(i2) '}$']);
%     zlabel(['$\theta_{' num2str(i3) '}$']);
%     legend('Estimates','Mean estimate','True value','$95\%$ confidence region','$95\%$ covariance ellipsoid')
%     legend('Location','northeast')
% %%    
% Theta_hor = [Theta_all(i1,1:10);Theta_all(i2,1:10);Theta_all(i3,1:10)];
%     Theta_hor = Theta_hor';
%     covar(:,1) = COV10([i1,i2,i3],i1);
%     covar(:,2) = COV10([i1,i2,i3],i2);
%     covar(:,3) = COV10([i1,i2,i3],i3);
%     [X_el,Y_el,Z_el] = build_ellipsoid(Theta_hor,covar);
%     [X_el1,Y_el1,Z_el1] = build_ellipsoid_from_data(Theta_hor);
%     fig('10',visFlag); 
%     plot3(Theta_hor(:,1),Theta_hor(:,2),Theta_hor(:,3),'o','Color',[80/255, 80/255, 80/255]); hold on;
%     plot3(Theta_mean_10(i1),Theta_mean_10(i2),Theta_mean_10(i3),'*r','LineWidth',5); hold on;
%     plot3(Theta_model(i1),Theta_model(i2),Theta_model(i3),'*k','LineWidth',5); hold on;
%     h = surf(X_el,Y_el,Z_el); alpha 0.8    
%     set(h, 'facecolor',[243/255, 222/255, 187/255]);
%     set(h, 'edgecolor',[237/255, 177/255, 32/255]);
%     hold on;
%     h1 = surf(X_el1,Y_el1,Z_el1); alpha 0.5    
%     set(h1, 'facecolor',[187/255, 222/255, 243/255]);
%     set(h1, 'edgecolor',[32/255, 177/255, 237/255]);
%     grid on;
%     view(45, 25);
% %     xlim([-100,150]); ylim([-100,150]); zlim([-100,150]);
%     xlabel(['$\theta_{' num2str(i1) '}$']);
%     ylabel(['$\theta_{' num2str(i2) '}$']);
%     zlabel(['$\theta_{' num2str(i3) '}$']);
%     legend('Estimates','Mean estimate','True value','$95\%$ confidence region','$95\%$ covariance ellipsoid')
%     legend('Location','northeast')
%% Pair-wise plots - uncomment if needed nore detail
% for iPair=1:2:ll-1
%     cc(:,1) = COV100([iPair, iPair+1],iPair);
%     cc(:,2) = COV100([iPair, iPair+1],iPair+1);
%     co{iPair} = cc;
%     mea{iPair} = [mean(Theta_all(iPair,:)),mean(Theta_all(iPair+1,:))];
%     ell_points_theta{iPair} = build_ellips(mea{iPair},co{iPair});
% end
% 
% for iPair=1:2:ll-1
%     fig(['Pair ',num2str(iPair)],visFlag); 
%     plot(Theta_all(iPair,:),Theta_all(iPair+1,:),'*k'); hold on;
%     plot(Theta_model(iPair),Theta_model(iPair+1),'*b','Linewidth',3); hold on;
%     plot(mea{iPair}(1),mea{iPair}(2),'*r','Linewidth',3); hold on;
%     plot(ell_points_theta{iPair}(:,1),ell_points_theta{iPair}(:,2),'r');
%     legend('MLEs','Mean MLE',)
%     xlabel(['$\Theta_{' num2str(iPair) '}$']);
%     ylabel(['$\Theta_{' num2str(iPair+1) '}$']);
% end
