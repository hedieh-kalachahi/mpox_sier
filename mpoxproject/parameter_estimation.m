% modeling for 120 days
t = [0:1:120];
%par=[S0,E0,I0,R0,beta,lambda,gamma], S0=339,996,563-10, E0=5, I0=5,R0=0,
%beta=0.00006, lambda=0.2, gamma=0.83
par=[339996553, 5,5,0,0.00006,0.2,0.83];
%% Parameters to Fit

% beta,lambda, gamma and p

widths = [2.0e6 2.0e4; 0.01 1; 1 5; 1e-4 1];

%Prior for beta,lambda, gamma are the following lognormal 
%distributions and for p is uniform distribution

prior1 = @(beta) lognpdf(log(beta),sigma=0.7);
prior2 = @(lambda) lognpdf(log(lambda),sigma=0.7);
prior3 = @(gamma) lognpdf(log(gamma),sigma=0.7);
prior4 = @(p) unifpdf(p, widths(4,1), widths(4,2));
par_names = ["beta", "lambda", "gamma", "p"];
%% Parameter Fit Setup
logl = @(theta) loglike_mpox_model(theta, data);
logp = @(theta) log(prior1(theta(1))) + log(prior2(theta(2))) + log(prior3(theta(3))) + log(prior4(theta(4)));

theta0 = zeros(4, 8); 
%Generate the vector of initial conditions

k = 0;
j = 1;

while(j < (8+1))

while(k == 0 || isnan(k) || isinf(k))
    
    initial = [widths(1,1)+(widths(1,2)-widths(1,1))*rand widths(2,1)+(widths(2,2)-widths(2,1))*rand widths(3,1)+(widths(3,2)-widths(3,1))*rand widths(4,1)+(widths(4,2)-widths(4,1))*rand];
    k = logl(initial) + logp(initial);
    
end

theta0(1,j) = initial(1);
theta0(2,j) = initial(2);
theta0(3,j) = initial(3);
theta0(4,j) = initial(4);

j = j + 1;
k = 0;

display(j);

end

save theta0.mat theta0
logfuns = {@(theta)logp(theta) @(theta)logl(theta)};
%% Affine Fitting Algorithm

%'gwmcmc' runs the MCMC sampler with affine invariance
[models, logP] = gwmcmc(theta0, logfuns, 1000000);

%'models' holds the proposed guesses for the parameter values
save Fitting_Example_models.mat models

%'logP' holds the corresponding value of the log prior and the value of the
%log likelihood for each guess of the parameter values
save Fitting_Example_logP.mat logP
%% Assess results

%Let 'l' hold the values of the log prior and log likelihood
l = logP;

%Set each of the eight chains to have a burn-in of 2500 samples
l(:,:,1:2500)=[];

%Extract the values of the log prior
lprior = l(1,:,:);

%Collapse the samples of all chains into a single larger sample
lprior = lprior(:,:);

%Extract the values of the log likelihood
llik = l(2,:,:);

%Collapse the samples of all chains into a single larger sample
llik = llik(:,:);

%Add the log prior and log likelihood together to get the log unnormalized 
%posterior samples
l_un_post_samples = lprior + llik;

%Let 'm' hold the values of the proposed guesses for the parameter values
m = models;

%Set each of the eight chains to have a burn-in of 2500 samples
m(:,:,1:2500)=[];

%Collapse the samples of all chains into a single larger sample
m = m(:,:);

theta_samples = m;

%'k1' contains all the proposed values for parameter 'x0'
k1 = m(1,:);

%'k2' contains all the proposed values for parameter 'y0'
k2 = m(2,:);

%'k3' contains all the proposed values for parameter 'log10_beta'
k3 = m(3,:);

%'k4' contains all the proposed values for parameter 'p'
k4 = m(4,:);

%parameter values at the maximum log unnormalized posterior value
[M, I] = max(l_un_post_samples);
%% Posterior Predictive distribution Loop

%The following 'for loop' estimates the posterior predictive distribution
%for the model solution x

x_predict = zeros(length(data),length(l_un_post_samples));
y_predict = zeros(length(data),length(l_un_post_samples));

discrepancy = zeros(1,length(l_un_post_samples));

discrepancy_pred = zeros(1,length(l_un_post_samples));

ind_D_pred_exceeds = zeros(1,length(l_un_post_samples));

h = waitbar(0,'Initialize...');
for i = 1:length(l_un_post_samples)
    
    theta = theta_samples(:,i);
    
    out = mpox_model(theta(1), theta(2), theta(3));
    
    discrepancy(1,i) = sum(((data(1,:) - out(1,:)).^2)./(out(1,:)/(theta(4)))) + sum(((data(2,:) - out(2,:)).^2)./(out(2,:)/(theta(4))));
    
    x_predict(:,i) = multnbinrnd(out(1,:)', theta(4), length(data(1,:)));
    y_predict(:,i) = multnbinrnd(out(2,:)', theta(4), length(data(2,:)));
    
    discrepancy_pred(1,i) = sum((((x_predict(:,i))' - out(1,:)).^2)./(out(1,:)/(theta(4)))) + sum((((y_predict(:,i))' - out(2,:)).^2)./(out(2,:)/(theta(4))));
    
    if(discrepancy_pred(1,i) > discrepancy(1,i))
        
        ind_D_pred_exceeds(1,i) = 1;
        
    else
        
        ind_D_pred_exceeds(1,i) = 0;
        
    end
    
   waitbar(i/length(l_un_post_samples),h,sprintf('%d%%',(i/length(l_un_post_samples))*100))
   
end
close(h)

save x_predict_mpox.mat x_predict
save y_predict_mpox.mat y_predict
save discrepancy_mpox.mat discrepancy
save discrepancy_pred_mpox.mat discrepancy_pred
save ind_D_pred_exceeds_mpox.mat ind_D_pred_exceeds

