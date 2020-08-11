function[X_out,P_out,Mu_out,X_cond,P_cond] = IMM_output(T,n_models,x_merged,P_merged,x_s,P_s,mu_s);
    for t=1:T
        X_out(:,t)   = x_merged{t};    % merged state
        P_out{t}(:,:) = P_merged{t};  % merged covariance
        Mu_out(:,t)  = mu_s{t}';     % mode probability
        for j=1:n_models
            % output to be used in the maximisation
            X_cond{j}(:,t)    = x_s{t,j};   % mode-matched posterior estimate
            P_cond{j,t}(:,:)  = P_s{t,j};   % mode-matched posterior covariance
        end
    end % for IMM assignment of values (t)