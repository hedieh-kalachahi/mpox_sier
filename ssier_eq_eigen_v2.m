% Define the system of differential equations
syms S E I R S_e E_e I_e

N = 15000;
b = 0.00003333;
beta = 0.00006;
beta_e = beta-beta*0.17;
alpha = 0.05;
sigma = 0.2;
gamma = 0.83;

dSdt = b*N-b*S-alpha*S-beta*S*I;
dS_edt = alpha*S-b*S_e-beta_e*S_e*I;
dEdt = beta*S*I-b*E-alpha*E-sigma*E;
dE_edt = beta_e*S_e*I+alpha*E-b*E_e-sigma*E_e;
dIdt = sigma*E-b*I-alpha*I-gamma*I;
dI_edt = sigma*E_e+alpha*I-b*I_e-gamma*I_e;
dRdt = gamma*I+gamma*I_e-b*R;

% Find the critical points by setting all derivatives equal to zero
critical_points = solve(dSdt == 0, dS_edt == 0, dEdt == 0, dE_edt == 0, dIdt == 0, dI_edt == 0, dRdt == 0, 'Real', true);

% Convert the critical points from symbolic to numeric values
cp = double([critical_points.S, critical_points.E, critical_points.I, critical_points.R, critical_points.S_e, critical_points.E_e, critical_points.I_e])

% Calculate the Jacobian matrix at each critical point, find its eigenvalues, and determine stability
for i = 1:size(cp,1)
    J = jacobian([dSdt, dS_edt, dEdt, dE_edt, dIdt, dI_edt, dRdt], [S, E, I, R, S_e, E_e, I_e]);
    J_num = double(subs(J, [S, E, I, R, S_e, E_e, I_e], cp(i,:)));
    disp(['Jacobian matrix at critical point ', num2str(i), ':']);
    disp(J_num);
    eig_J = eig(J_num);
    disp(['Eigenvalues of Jacobian matrix at critical point ', num2str(i), ':']);
    disp(eig_J);
    if all(real(eig_J) < 0)
        disp('All eigenvalues have negative real part, the critical point is stable');
    elseif all(real(eig_J) <= 0) && any(real(eig_J) == 0)
        disp('All eigenvalues have nonpositive real part, but at least one eigenvalue is zero, further analysis is needed');
    else
        disp('At least one eigenvalue has positive real part, the critical point is unstable');
    end
    disp('---------------------------------------');
end

