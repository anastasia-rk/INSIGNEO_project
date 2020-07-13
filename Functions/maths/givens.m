function [Rotation,M] =  givens(M)
    for i=1:size(M,2)
        for j =size(M,1):-1:i+1
            [G,M] = rotate(M,i,j);
        end
    end
    Rotation = G;
end

function[G,M] = rotate(M,i,j)
% M - matrix that is subjected to Givens rotations
% i - column of interest
% j - row of interest
G = eye(size(M,1));
if (M(j,i) ~= 0) 
    v = M(j-1:j,i);
    G_v = planerot(v);
    G(j-1,j-1) = G_v(1,1);
    G(j-1,j) = G_v(1,2);
    G(j,j-1) = G_v(2,1);
    G(j,j) = G_v(2,2);
    M = G*M;
end
end
