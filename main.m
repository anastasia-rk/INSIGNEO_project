local_init;
%% Select the pattern
nTracks = 100;                                                              % number of tracks to generate in each map
pattern = questdlg('Starting position distribution', ...
    'Select pattern',...
	'Uniform','Line','Point','');
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
nx = 4; ny = 4;
switch basis_type
    case 'gaussian'
        [knots,sigma] = setup_gaussian_support(grid_limits,nx,ny);
        % equal sigmas for isotopic basis functions
        sigma = 80;
        Z = [sigma^2 0; 0 sigma^2];   
        ll = size(knots,2);
    case 'bspline'
        [knots] = setup_spline_support(grid_limits,nx,ny);
        Z = 0;
        ll = size(knots,2)/2;
end
%% Initialise field estimate and model parameters
grid_limits1 = grid_limits;
Theta = initiate_field(0,60,grid_limits1,knots);                            %fit splines to a flat surface
Theta_old = Theta;

 Theta_model = ones(ll,1);
 Theta = ones(ll,1);
 Theta_model(1:4,1) = 10;
 Theta_model(5:8,1) = 100;
 Theta_model(9:12,1) = 190;
 Theta_model(13:16,1) = 290;
 
figure; 
colormap(my_map);
done = plot_field(Theta_model,Z,knots,grid_limits,'bspline');
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
H = [I O];
% Disturbance matrix
G_cv = [T^2/2*I; T*I];
G_rw = [T^2/2*I; T*I]; 
%% Covariances for the Kalman filters
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
% for the maximisation step
Sig_w{1} = G{1}*(inv(G{1}'*G{1}))'*Q_cv_inv*inv(G{1}'*G{1})*G{1}';
Sig_w{2} = G{2}*(inv(G{2}'*G{2}))'*Q_rw_inv*inv(G{2}'*G{2})*G{2}';
%% Monte Carlo loop
models = [1:50];                                                            % Monte-Carlo simulations
for iModel = models
clear X Y Mode_model
load([folderName,'/simulated_tracks_' num2str(iModel)]);
Tracks = [1:100]; 
% Tracks = 1; % for now
converged_theta = false;
converged_phi   = false;
converged_r     = false;
iter_max  = 10;
iter      = 0;
Theta = zeros(ll,1); 
for k=Tracks
    Z_temp{k} = zeros(4,length(Y{k}));
    Z_temp{k}(1:2,:) = Y{k};
end
%% EM loop
while  (iter < iter_max)
%% Expectation step
fprintf('E-step... \n')
clear x_merged p_merged mu_s mu x_m x_s P_m P_s
iter = iter + 1;  
for k=Tracks    
T = length(Y{k});
for t=1:T
    y{t} = Y{k}(:,t);
end
%% recurtion cycle for IMM filter
mu{1} = mu_0;
for j=1:n_models
    x_m{1,j} = zeros(x_len,1);
    x_m{1,j}(1:2,:) = y{1};
    P{1,j} = eye(x_len);
end
for t=2:T
    % Mode conditioned calculations
   for j=1:n_models
       % Calculate mixing probabilities
       c(j) = mu{t-1}*p_tr(:,j); % normalising constant
       for i=1:n_models 
           mu_mix(i,j) = p_tr(i,j)*mu{t-1}(i)/c(j); % probabilities
       end
       % Mixing initial conditions for filtering
       v = zeros(x_len,1);
       for i=1:n_models
           v = v +  x_m{t-1,i}*mu_mix(i,j);
       end
       pv = zeros(x_len);
       for i=1:n_models
           pv = pv +  mu_mix(i,j)*(P{t-1,i} + (x_m{t-1,i} - v)*(x_m{t-1,i} - v)');
       end
       x_0{j} = v;
       P_0{j} = pv;
       % Model-matched filtering
       % Linear
       clear beta;
%        beta = field_gradient(Z_temp{k}(1:2,t-1),Z,knots,basis_type);
%        u{t-1,j} = mu_field*beta*Theta;
%        [x_m{t,j},P{t,j},lik{t,j}] = kf(y{t},x_0{j},u{t-1,j},P_0{j},F{j},B{j},C,Q{j},R);
       % Nonlinear
         [x_m{t,j},P{t,j},lik{t,j}] = ekf(y{t},x_0{j},P_0{j},Q{j},R,F{j},B{j},C,Theta,Z,knots,basis_type);
%         [x_pr,x_m{t,j},P_pr,P{t,j},lik{t,j}] = ukf(y{t},x_0{j},P_0{j},Q{j},R,F{j},B{j},C,Theta,Z,knots,basis_type);            
       % Prior to probabilities update
       m(1,j) = lik{t,j}*c(j);
       lik_plot(j,t) = lik{t,j};
       x_cond_plot{j}(:,t) = x_m{t,j};
   end
   % Mode probabilites 
   c_all = sum(m);
   mu{t} = m./c_all;
   % Merging states (optional)
   v = zeros(x_len,1);
   for j=1:n_models
       v = v + mu{t}(1,j)*x_m{t,j};
   end
   x_merged{t} = v;
   % Merging covariances (optional)
   pv = zeros(x_len);
   for j=1:n_models
       pv = pv + mu{t}(1,j)*(P{t,j}+ (x_m{t,j} - v)*(x_m{t,j} - v)');
   end
   P_merged{t} = pv;
   x_plot(:,t) = v;
   mu_plot(:,t) = mu{t}';
   [max_mode,i] = max(mu{t});
   mode_probable{k}(t) = i;
   X_filt{k}(:,t) = x_merged{t};    % merged state
   Mu_filt{k}(:,t)  = mu{t}'; 
end
   X_out{k}(:,t) = x_merged{T};    % merged state
   P_out{k,T}(:,:) = P_merged{T};  % merged covariance
   P_for_resampling{k}(T,:,:) = P_merged{T};  % merged covariance
   for j =1:n_models
        X_cond{k,j}(:,T) = x_m{T,j};     % mode-matched posterior estimate
        P_cond{k,j,T}(:,:) = P{T,j};   % mode-matched posterior covariance
   end
%% Recurtion cycle for IMM RTS smoother
% Initialise smoother
mu_s{T} = mu{T};
for j=1:n_models
    x_s{T,j} = x_m{T,j};
    P_s{T,j} = P{T,j};
end
   Mu_out{k}(:,T)  = mu_s{T}';     % mode probability

for t=T-1:-1:1 
% Backward transition kernel
    for i=1:n_models
       % Calculate mixing probabilities
       e(i) = mu{t}*p_tr(:,i); % normalising constant
       for j=1:n_models 
           b_tr(i,j) = p_tr(j,i)*mu{t}(j)/e(i); % probabilities
       end
    end
    % Mode conditioned calculations
   for j=1:n_models
       % Calculate mixing probabilities
       d(j) = mu_s{t+1}*b_tr(:,j); % normalising constant
       for i=1:n_models 
           mu_mix(j,i) = b_tr(i,j)*mu_s{t+1}(i)/d(j); % probabilities
       end
       mu_cond{k,t} = mu_mix;
       % Mixing initial conditions for filtering
       v = zeros(x_len,1);
       for i=1:n_models
           v = v +  x_s{t+1,i}*mu_mix(j,i);
       end
       pv = zeros(x_len);
       for i=1:n_models
           pv = pv +  mu_mix(j,i)*(P_s{t+1,i} + (x_s{t+1,i} - v)*(x_s{t+1,i} - v)');
       end
       xs_0{j} = v;
       Ps_0{j} = pv;
       % Model-matched smoothing
       % Linear
%        beta = field_gradient(Z_temp{k}(1:2,t+1),Z,knots,basis_type);
%        u{t,j} = mu_field*beta*Theta;
%        [x_s{t,j},P_s{t,j},lik_s{t,j}] =       rts(xs_0{j},x_m{t,j},u{t,j},Ps_0{j},P{t,j},F{j},B{j},Q{j});
        [x_s{t,j},P_s{t,j}] = erts(xs_0{j},x_m{t,j},Ps_0{j},P{t,j},F{j},B{j},Q{j},Theta,Z,knots,basis_type); % ,lik_s{t,j}
%          [x_s{t,j},P_s{t,j},lik_smm] = urts(xs_0{j},x_m{t,j},Ps_0{j},P{t,j},F{j},B{j},G{j},Q{j},Theta,Z,knots,basis_type);
       
       % Priot to probabilities update
       x_sm_cond_plot{j}(:,t) = x_s{t,j};
       % output to be used in the maximisation
       X_cond{k,j}(:,t) = x_s{t,j};     % mode-matched posterior estimate
       P_cond{k,j,t}(:,:) = P_s{t,j};   % mode-matched posterior covariance
   end
   % Compute smoothed likelihoods
   for j=1:n_models
       lik_s{t,j} = 0;
       for i=1:n_models
           n_s{t+1,i} = sm_lik(x_s{t+1,i},x_m{t,j},P{t,j},F{j},B{j},Q{j},Theta,Z,knots,basis_type);
           lik_s{t,j} = lik_s{t,j} + n_s{t+1,i}*p_tr(i,j);
           lik_s_plot(j,t) = lik_s{t,j};
       end
       m_s(1,j) = lik_s{t,j}*mu{t}(1,j);
   end
% Mode probabilites 
   d_all = sum(m_s);
   mu_s{t} = m_s./d_all;
% Merging states
   v = zeros(x_len,1);
   for j=1:n_models
       v = v + mu_s{t}(1,j)*x_s{t,j};
   end
   x_merged{t} = v;
% Merging covariances
   pv = zeros(x_len);
   for j=1:n_models
       pv = pv + mu_s{t}(1,j)*(P{t,j}+ (x_s{t,j} - v)*(x_s{t,j} - v)');
   end
   P_merged{t} = pv;
   x_sm_plot(:,t) = v;
   mu_plot(:,t) = mu_s{t}';
   [max_mode,i] = max(mu_s{t});
   mode_probable{k}(t) = i;
   
% Output of the IMM algorithm
   Lik_out{k}(:,t) = lik_s_plot(:,t);
   X_out{k}(:,t) = x_merged{t};    % merged state
   P_out{k,t}(:,:) = P_merged{t};  % merged covariance
   P_for_resampling{k}(t,:,:) = P_merged{t};  % merged covariance
   Mu_out{k}(:,t)  = mu_s{t}';     % mode probability
end % for smoother (t)
Z_temp{k} = X_out{k};
end % for track (k)

%% Maximization step [2]
fprintf('M-step... \n')
%% Transition probability matrix 
mu_joint = zeros(n_models,n_models);
for j=1:n_models
    denom(j) = 0;
    for l=1:n_models
        for k=Tracks
            [nn,T] = size(X_out{k});
            for t=1:T-2
                mu_joint(j,l) = mu_joint(j,l) + Mu_out{k}(j,t)*mu_cond{k,t+1}(j,l);
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
sum1 = 0; sum2 = 0; sum3 = 0;
denom4 = 0; sum4 = zeros(2,2); sum5 = zeros(2,2);
for k=Tracks
      clear dx x1 x2 y x
      y = Y{k};
      [nn,T] = size(X_out{k}); 
      x = X_out{k};
      % CV mode
      x1(:,:) =  X_cond{k,1};
      dx{1}(:,1:T-1) = x1(:,2:T) - F_cv*x(:,1:T-1);
      dx{1}(:,T) = dx{1}(:,T-1);
      % RW mode
      x2(:,:) =  X_cond{k,2};
      dx{2}(:,1:T-1) = x2(:,2:T) - F_rw*x(:,1:T-1);
      dx{2}(:,T) = dx{2}(:,T-1);
     for t=1:T-1
         gg(:,:) = field_gradient(Z_temp{k}(1:2,t),Z,knots,basis_type);
         for j=1:n_models
         bb = B{j}*mu_field*gg;
         sum1 = sum1 + Mu_out{k}(j,t)*bb'*Sig_w{j}*bb;
         sum2 = sum2 + Mu_out{k}(j,t)*bb'*Sig_w{j}*dx{j}(:,t);
         sum3 = sum3 + Mu_out{k}(j,t)*bb'*Sig_w{j}*bb;
         end % for modes (j) 
         sum5 = sum5 + y(:,t)*y(:,t)' - y(:,t)*x(:,t)'*C' - C*x(:,t)*y(:,t)' + C*(P_out{k,t} + x(:,t)*x(:,t)')*C';
     end % for times (t)
      sum4 = sum4 + (y - C*x)*(y - C*x)';
    denom4 = denom4 + T;
end % for tracks (k)
% R_old = R;
% R = sum5/denom4
R_est{iter} = sum5/denom4;
clear y     
% Parameter MLE
Theta_old = Theta;
Theta = pinv(sum1)*sum2;
FieldCoeff{iter} = Theta;
% Expected Fisher information
Expected_info = sum3;
Missing_info = (sum1*Theta - sum2)*(sum1*Theta_old - sum2)';
% temp(:,iter) = Theta;
%% Check convergence
fprintf('Checking convergence condition... \n')
[converged_theta]   = converge_param(Theta,Theta_old,iter);
[converged_phi]     = converge_param(vec(p_tr),vec(p_tr_old),iter);
[converged_r]       = converge_param(vec(R),vec(R_old),iter);
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
[done1,ZZ] = plot_heatmap(Theta_diff,Z,knots,grid_limits,basis_type);
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
[done2,ZZ] = plot_heatmap(Theta_mean,Z,knots,grid_limits,basis_type);
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
[done3, ZZ1] = plot_gradient(Theta_model,Z,knots,grid_limits,basis_type);
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
[done3, ZZ1] = plot_gradient(Theta_mean,Z,knots,grid_limits,basis_type);
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
done4 = plot_gradient(Theta,Z,knots,grid_limits,basis_type);
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
done5 = plot_gradient(Theta,Z,knots,grid_limits,basis_type);
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
