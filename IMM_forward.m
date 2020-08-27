function[x_m,P,mu,mode_probable] = IMM_forward(T,x_len,n_models,p_tr,mu_0,Y,Q,R,F,B,C,Theta,knots,order)
for t=1:T
    y{t} = Y(:,t);
end
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
         [x_m{t,j},P{t,j},lik{t,j}] = ekf(y{t},x_0{j},P_0{j},Q{j},R,F{j},B{j},C,Theta,knots,order);
%         [x_pr,x_m{t,j},P_pr,P{t,j},lik{t,j}] = ukf(y{t},x_0{j},P_0{j},Q{j},R,F{j},B{j},C,Theta,Z,knots,basis_type);            
       % Prior to probabilities update
       m(1,j) = lik{t,j}*c(j);
       lik_plot(j,t) = lik{t,j};
       x_cond_plot{j}(:,t) = x_m{t,j};
   end
   % Mode probabilites 
   c_all = sum(m);
   mu{t} = m./c_all;
      [max_mode,i] = max(mu{t});
   mode_probable(t) = i;
%%
%    % Merging states (optional)
%    v = zeros(x_len,1);
%    for j=1:n_models
%        v = v + mu{t}(1,j)*x_m{t,j};
%    end
%    x_merged{t} = v;
%    % Merging covariances (optional)
%    pv = zeros(x_len);
%    for j=1:n_models
%        pv = pv + mu{t}(1,j)*(P{t,j}+ (x_m{t,j} - v)*(x_m{t,j} - v)');
%    end
%    P_merged{t} = pv;
end % for filter (t)
end