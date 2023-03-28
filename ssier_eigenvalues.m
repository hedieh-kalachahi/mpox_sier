% Define the parameters of the modified SIER model
beta = 0.5;
sigma = 0.1;
gamma = 0.05;
N = 1000;
p = 0.1;

% Define the Jacobian matrix
J = [-beta*(I+E)/(N+E) - beta*S*(I+E)/(N+E)^2, -beta*S/(N+E), 0, 0, 0, 0;
     beta*(I+E)/(N+E) + beta*S*(I+E)/(N+E)^2, beta*S/(N+E) - sigma, 0, 0, 0, 0;
     0, sigma, -(1-p)*gamma, 0, -p*gamma, 0;
     0, 0, (1-p)*gamma, 0, 0, 0;
     0, 0, 0, 0, -sigma, p*gamma;
     0,];
% Compute the eigenvalues of the Jacobian matrix
lambda = eig(J);

% Print the eigenvalues
disp(lambda);