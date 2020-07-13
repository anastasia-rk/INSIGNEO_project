function[x_t1t1,P_t1t1,lik] = kf(y,x_tt,u,P_tt,A,B,C,Q,R)
% tt = (t \mid t)
% t1t = (t + 1 \mid t)
% t1t1 = (t + 1 \mid t + 1)
%% Estimation
x_t1t   = A*x_tt + B*u;                                                     % prior mean
P_t1t   = A*P_tt*A' + Q;                                                    % prior covariance
d_y     = y - C*x_t1t;                                                      % estimation error    
S       = C*P_t1t*C' + R;                                                   % error covariance
S_inv   = pinv(S);
K       = P_t1t*C'*S_inv;                                                   % Kalman gain
x_t1t1  = x_t1t + K*d_y;                                                    % postrerior mean          
P_t1t1  = (eye(length(x_tt)) - K*C)*P_t1t;                                  % posterior covariance
%% model likelihood 
S_new   = (abs(2*pi*S));
den     = sqrt(det(S_new));
num     = exp(-0.5*d_y'*S_inv*d_y);
lik     = num/den;