% Define the parameters of the SIER model
beta = 0.5;
sigma = 0.1;
gamma = 0.05;
N = 1000;

% Define the Jacobian matrix
J = [-beta/N, -beta/N, 0, 0;
     beta/N, beta/N-sigma, 0, 0;
     0, sigma, -gamma, 0;
     0, 0, gamma, 0];

% Compute the eigenvalues of the Jacobian matrix
lambda = eig(J);

% Print the eigenvalues
disp(lambda);

