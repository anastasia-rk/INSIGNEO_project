function[x_tT,P_tT] = rts(x_t1T,x_tt,P_t1T,P_tt,u,A,B,Q)
% tt = (t \mid t)
% t1t = (t + 1 \mid t)
% t1T = (t + 1 \mid T)
% tT = (t \mid T)
%% Estimation
x_t1t   = A*x_tt + B*u;                                                  % prior estimate from Kalman filter
P_t1t   = A*P_tt*A' + Q;                                                 % prior covariance from KF
Gain    = P_tt*A'*pinv(P_t1t);                                                 % smoother gain
x_tT    = x_tt + Gain*(x_t1T - x_t1t);
P_tT    = P_tt + Gain*(P_t1T - P_t1t)*Gain';
%% Model likelihood
d_x     = x_t1T - x_t1t;
S_new   = (abs(2*pi*P_t1t));
den     = sqrt(det(S_new));
num     = exp(-0.5*d_x'*pinv(P_t1t)*d_x);
lik     = num/den;