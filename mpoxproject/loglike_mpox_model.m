function [loglike] = loglike_mpox_model(theta, data)


beta = theta(1);
lambda = theta(2);
gamma = theta(3);
p = theta(4);
par=[339996553, 5,5,0,0.00006,0.2,0.83];
time=120
m_out = mpox_model(par, time);

check_size = size(m_out);


if(check_size(1) == length([0:1:120]))

loglike_1 = sum(log(nbinpdf(data(1,:), (m_out(1,:)*p)/(1 - p), p)));

loglike = loglike_1;

else

    loglike = -Inf;

end

end