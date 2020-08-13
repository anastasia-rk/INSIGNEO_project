function[x_tT,P_tT,lik] = erts(x_t1T,x_tt,P_t1T,P_tt,A,B,Q,theta,knots,order)
% tt = (t \mid t)
% t1t = (t + 1 \mid t)
% t1T = (t + 1 \mid T)
% tT = (t \mid T)
%% Taylor expansion
dx = 0.001;
for i=1:length(x_tt)
    xx      = x_tt;
    xx(i)   = xx(i) + dx; 
    F(:,i) = (dynfun (xx,A,B,theta,knots,order) - dynfun (x_tt,A,B,theta,knots,order))./dx;
end
%% Suboptimal estimation
x_t1t   = dynfun (x_tt,A,B,theta,knots,order);                               % prior estimate from Kalman filter
P_t1t   = F*P_tt*F' + Q;                                                    % prior covariance from KF
Gain    = P_tt*F'*pinv(P_t1t);                                              % smoother gain
x_tT    = x_tt + Gain*(x_t1T - x_t1t);
P_tT    = P_tt + Gain*(P_t1T - P_t1t)*Gain';   
P_tT    = (P_tT + P_tT.')/2;
%% Model likelihood
d_x     = x_t1T - x_t1t;
S_new   = (abs(2*pi*P_t1t));
den     = sqrt(det(S_new));
num     = exp(-0.5*d_x'*pinv(P_t1t)*d_x);
lik     = num/den;