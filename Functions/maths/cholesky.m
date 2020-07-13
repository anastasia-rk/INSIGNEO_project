function[S] = cholesky(P)
n = size(P,2);
S = zeros(n);
for i=1:n
    a = 0;
    for j=1:i-1
        a = a + S(i,j)^2;
    end
    S(i,i) = sqrt(P(i,i) - a);
    for j=i:n
        b = 0;
        for k=1:i-1
            b = b + S(j,k)*S(i,k);
        end
        S(j,i) = (P(j,i) - b)/S(i,i);
    end
end
end