function[mixture] = gm_elements(x_s,P_s,mu_s)
%%
eps = 10e-8;
[T,N] = size(x_s);
for t =1:T
    mu_mixture = mu_s{t};
    for j=1:N
        mean_mixute(:,j)    = x_s{t,j};  
        P_mixture(:,:,j)    = (P_s{t,j} + P_s{t,j}.')/2; 
        %% mixture weights can't be zero - replace by an eps weight
        if mu_mixture(j) < eps 
            mu_mixture(j) = eps;
        end
    end
    % matlab renormalises mixture weights in case they don't sum to one
    mixture{t} = gmdistribution(mean_mixute',P_mixture,mu_mixture);
end
