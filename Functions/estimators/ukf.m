function[x_t1t1,P_t1t1,lik] = ukf(y,x_tt,u,P_tt,dfun,mfun,Q,R)
% tt = (t \mid t)
% t1t = (t + 1 \mid t)
% t1t1 = (t + 1 \mid t + 1)
%% Preliminaries
L       = length(x_tt);                                                     % state vector dimension
M       = length(y);                                                        % measurement vector dimension
alpha   = 0.01;                                                             % scaling constant - from literature (set to larger value for wider distribution spread)
kappa   = 0;                                                                % from literature
beta    = 2;                                                                % optimal in case distribution of x is Gaussian
lambda  = alpha^2*(L + kappa) - L;                                          % scaling factor for sigma points
c       = L + lambda;                                                       % scaling factor for covariance
x_t1t = zeros(L,1);
y_t1t = zeros(M,1);
P_t1t = zeros(L,L);
P_yy    = zeros(M,M);
P_xy    = zeros(L,M);
%% Sigma ponts
[LM,DM,E,pneg]  = mchol(P_tt);
S               = sqrt(c)*LM*sqrt(DM);
X{1}            = x_tt;
xxx_t           = x_tt(:,ones(1,numel(x_tt)));
X3              = [x xxx_t+S xxx_t-S]; 
for i = 2:2*L+1
    X{i} = X3(:,i);
end
% Weights
W_s(1)          = lambda/c;
W_c(1)          = lambda/c + (1 - alpha^2 + beta);
W_s (2:2*L+1)   = 1/(2*c);                                                  % for states
W_c (2:2*L+1)   = 1/(2*c);                                                  % for covariances
%% Unscented transform
for i=1:2*L+1
    Xx{i}   = dfun(X{i},u); % sigma state propagation
    x_t1t   = x_t1t + W_s(i)*Xx{i};                                         % prior estimate
    Yy{i}   = mfun(Xx{i},u);                                                % sigma measuremnent
    y_t1t   = y_t1t + W_s(i)*Yy{i};                                         % measurement of prior
end

for i=1:2*L+1
    P_t1t   = P_t1t + W_c(i)*(Xx{i} - x_t1t)*(Xx{i} - x_t1t)' + Q;          % prior covariance
    P_yy    = P_yy + W_c(i)*(Yy{i} - y_t1t)*(Yy{i} - y_t1t)' + R;           % residual covariance
    P_xy    = P_xy + W_c(i)*(Xx{i} - x_t1t)*(Yy{i} - y_t1t)';               % cross-covariance
end
%% Correction
K       = P_xy*pinv(P_yy);                                                  % Kalman gain
x_t1t1  = x_t1T + K*(y - y_t1t);                                            % state correction
P_t1t1  = P_t1t - K*P_yy*K';                                                % covariance correction
%% Model likelihood
S_new   = (abs(2*pi*S));
den     = sqrt(det(S_new));
num     = exp(-0.5*d_y'*S_inv*d_y);
lik     = num/den;