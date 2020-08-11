function[X_resample,P_resample,nResampled,indeces] = resample_tracks(X_smoothed,P_smoothed,X_mode,mode,nTracks)
X_resample = [];
P_resample = [];
indeces = [];
nResampled = nTracks;
jj = 0;
for j = 1:nTracks
    kk = 0; % for storing indeces
    indeces{j} = 0;
    clear index
    [n,m] = size(X_smoothed{j});
    index = find(X_mode{j} == mode);
    if length(index) < 2 % none/only one time point in a sensing segment
        nResampled = nResampled - 1;
    else % several time points in a track
        nn = 1; % number of segments in a given track
        i_b(nn) = 1;
        for i = 2:length(index) % 
            if (index(i) - index(i-1) > 1)
                nResampled = nResampled + 1; 
                nn = nn + 1;
                i_e(nn-1) = i-1;
                i_b(nn) = i;
            end
        end
        i_e(nn) = length(index);
        for k = 1:nn
        if (i_e(k) - i_b(k)) > 0
        jj = jj + 1;
        kk = kk + 1;
        indeces{j}(kk) = jj;
        X_resample{jj} = zeros(length(index(i_b(k):i_e(k))),n);
        P_resample{jj} = zeros(length(index(i_b(k):i_e(k))),4,4);
        X_resample{jj} = X_smoothed{j}(1:n,index(i_b(k):i_e(k)))';
        P_resample{jj} = P_smoothed{j}(index(i_b(k):i_e(k)),:,:);
    else
        nResampled = nResampled - 1;
    end
    end
    end    
end
end
