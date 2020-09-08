function[x_tilde,denoms] = univmultivar(X_out,P_out,L)
a_x = X_out'-3*sqrt(diag(P_out))';
b_x = X_out'+3*sqrt(diag(P_out))';
x_len = length(a_x);
for i=1:L
    x_tilde(i,:) = unifrnd(a_x,b_x,1,x_len);
    denoms(i) = prod(unifpdf(x_tilde(i,:),a_x,b_x));
end