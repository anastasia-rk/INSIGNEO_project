for i=1:3
    for j=1:120
        pd_check(j,i) = isPositiveDefinite(P_s{j,i});
    end
end



%%
function x=isPositiveDefinite(A)
%Function to check whether a given matrix A is positive definite
%Author Mathuranathan for https://www.gaussianwaves.com
%Returns x=1, if the input matrix is positive definite
%Returns x=0, if the input matrix is not positive definite
%Throws error if the input matrix is not symmetric
    %Test for positive definiteness
    x=1; %Flag to check for positiveness
    [~,x] = chol(A);
end