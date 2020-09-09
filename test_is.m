clc; clear; close all
%%
mu = [0; 2];
sigma(1,1,1) = 0.2;
sigma(1,1,2) = 0.2;
p = [0.1 0.9];
gm = gmdistribution(mu,sigma,p);
mu1 = p*mu;
sigma1 = 0;
for i=1:2
   sigma1 = sigma1 + p(i)*(sigma(:,:,i) + (mu(i) - mu1)*(mu(i) - mu1)');
end

figure;
fplot(@(x)reshape(pdf(gm,x(:)),size(x)),[-10 10]);
hold on;
x = [-10:.1:10];
y = normpdf(x,mu1,sigma1);
plot(x,y);

%%
% mu1 = mu1 + 0.3
L = 300;
 x_tilde = mvnrnd(mu1,sigma1,L);
 x_tilde_gm = random(gm,L);
 numers_gm = pdf(gm,x_tilde_gm);
 numers = pdf(gm,x_tilde);
 denoms =  mvnpdf(x_tilde,mu1,sigma1);
 weights = numers./denoms;
 figure;
 subplot(2,1,1)
 scatter(x_tilde,numers); hold on;
 scatter(x_tilde,denoms); hold on;
 subplot(2,1,2)
 scatter(x_tilde,weights); hold on;
 func = @(x) x^2 + x^3;
 
 sum1 = 0;
 sum2 = 0;
 for i=1:L
     sum1 = sum1 + func(x_tilde_gm(i))*numers_gm(i);
     sum2 = sum1 + func(x_tilde(i))*weights(i)*denoms(i);
 end
 
 exp_mc = sum1/L
 exp_is = sum2/L


