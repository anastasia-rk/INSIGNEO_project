function[x_s,P_s,mu_s,x_merged,P_merged,mu_cond,mode_probable] = IMM_backward(T,x_len,n_models,p_tr,mu,x_m,P,F,B,Q,Theta,knots,order)
% Initialise smoother
mu_s{T} = mu{T};
for j=1:n_models
    x_s{T,j} = x_m{T,j};
    P_s{T,j} = P{T,j};
end
%% Assign final merged value at t=T
   % Merging states (optional)
   v = zeros(x_len,1);
   for j=1:n_models
       v = v + mu{T}(1,j)*x_m{T,j};
   end
   x_merged{T} = v;
   % Merging covariances (optional)
   pv = zeros(x_len);
   for j=1:n_models
       pv = pv + mu{T}(1,j)*(P{T,j}+ (x_m{T,j} - v)*(x_m{T,j} - v)');
   end
   P_merged{T} = pv;
%% Main loop
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
       mu_cond{t} = mu_mix;
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
        [x_s{t,j},P_s{t,j},~] = erts(xs_0{j},x_m{t,j},Ps_0{j},P{t,j},F{j},B{j},Q{j},Theta,knots,order); % ,lik_s{t,j}
%          [x_s{t,j},P_s{t,j},lik_smm] = urts(xs_0{j},x_m{t,j},Ps_0{j},P{t,j},F{j},B{j},G{j},Q{j},Theta,Z,knots,basis_type);
       
       % Priot to probabilities update
%        x_sm_cond_plot{j}(:,t) = x_s{t,j};

   end
   % Compute smoothed likelihoods
   for j=1:n_models
       lik_s{t,j} = 0;
       for i=1:n_models
           n_s{t+1,i} = sm_lik(x_s{t+1,i},x_m{t,j},P{t,j},F{j},B{j},Q{j},Theta,knots,order);
           lik_s{t,j} = lik_s{t,j} + n_s{t+1,i}*p_tr(i,j);
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
%    x_sm_plot(:,t) = v;
%    mu_plot(:,t) = mu_s{t}';
   [max_mode,i] = max(mu_s{t});
   mode_probable(t) = i;
end % for smoother (t)

end