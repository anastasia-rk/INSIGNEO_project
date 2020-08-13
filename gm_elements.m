function[mixture] = gm_elements(x_s,P_s,mu_s)
[T,N] = size(x_s);
for t =1:T
    for j=1:N
        mean_mixute(:,j)    = x_s{t,j};  
        P_mixture(:,:,j)    = (P_s{t,j} + P_s{t,j}.')/2; 
    end
    mixture{t} = gmdistribution(mean_mixute',P_mixture,mu_s{t});
end
