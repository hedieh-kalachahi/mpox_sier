function [out] = mpox_model(par,time)
% S=X1, E=X2, I=X3, R=X4.
S0=par(1);
E0=par(2);
I0=par(3);
R0=par(4);
beta=par(5);
lambda=par(6);
gamma=par(7);
time=120
t = [0:1:time];

f = @(t,x) [-beta*x(1)*x(3); beta*x(1)*x(3)-lambda*x(2);lambda*x(2)-gamma*x(3);gamma*x(3)];

[~, xa] = ode45(f, t, [S0,E0,I0,R0]);

out(1,:) = (xa(:,1))';
out(2,:) = (xa(:,2))';
out(3,:) = (xa(:,3))';
out(4,:) = (xa(:,4))';
end