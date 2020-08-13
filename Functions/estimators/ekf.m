function[x_t1t1,P_t1t1,lik] = ekf(y,x_tt,P_tt,Q,R,A,B,C,theta,knots,order)
% tt = (t \mid t)
% t1t = (t + 1 \mid t)
% t1t1 = (t + 1 \mid t + 1)
%% Taylor expansion
dx = 0.001;
for i=1:length(x_tt)
    xx      = x_tt;
    xx(i)   = x_tt(i) + dx;   
    F(:,i) = (dynfun (xx,A,B,theta,knots,order) - dynfun (x_tt,A,B,theta,knots,order))./dx;
    H(:,i) = (measfun(xx,C) - measfun(x_tt,C))./dx;
end
%% Suboptimal estimation
x_t1t   = dynfun (x_tt,A,B,theta,knots,order);  % prior mean
P_t1t   = F*P_tt*F' + Q;                        % prior covariance
d_y     = y - measfun(x_tt,C);                  % estimation error    
S       = H*P_t1t*H' + R;                       % error covariance
S_inv   = pinv(S);
K       = P_t1t*H'*S_inv;                       % Kalman gain
x_t1t1  = x_t1t + K*d_y;                        % postrerior mean          
P_t1t1  = (eye(length(x_tt)) - K*H)*P_t1t;      % posterior covariance  
P_t1t1  = (P_t1t1 + P_t1t1.')/2;
%% Model likelihood
S_new   = (abs(2*pi*S));
den     = sqrt(det(S_new));
num     = exp(-0.5*d_y'*S_inv*d_y);
lik     = num/den;
