function [z, rel_err, val, overall_time] = PR_Cubic_Reg(z_true, A, b, init, tol, k_max)

% This program solves the phase retrieval problem by cubic regularization.
% The inner problem is handled as a QEP (quadratic eigenvalue problem).

if (nargin == 3)
    init = 0;
    tol = 10^(-4);
    k_max = 100;
elseif (nargin == 4)
    tol = 10^(-4);
    k_max = 100;
elseif (nargin == 5)
    k_max = 1000;
end

n = length(z_true);

% Initialization
Init_Start = tic;
if init==0
    z0 = 10*rand(n,1)-5 + 1i*(10*rand(n,1)-5);
else
    fprintf('Initializing...,');
    z0 = PR_Init(A,b);
    fprintf('OK! Finished. Time consumed: %f\n\n\n', toc(Init_Start));
end
Init_Time = toc(Init_Start);

% Some Constants
J = [eye(n) 1i*eye(n);
     eye(n) -1i*eye(n)]./2;
L2 = 20*max(sqrt(sum(abs(A).^2,1)));

% Initial Point of SIPM for subproblem
GEP_Init_Pt = randn(4*(2*n)+2,1);
GEP_Init_Pt = GEP_Init_Pt/norm(GEP_Init_Pt);

% Subproblem error tolerance mode: 
%     1 stands for constant tolerance; 
%     2 stands for diminising tolerance.
sub_tol_mode = 2;

start_time = tic;
k = 0;
rel_err = Inf;
z = z0;
val = PR_f(z,A,b);
err_vec = zeros(k_max,1);
while ((rel_err > tol) && (k <= k_max))
    Subprob_Start = tic;
    k = k + 1;
    sigma = L2;
    
    % Form gradient
    Grad_Start = tic;
    g = PR_Grad(z,A,b);
    real_g = 2*J'*g;
    Grad_Time = toc(Grad_Start);
    
    % Form Hessian
    Hess_Start = tic;
    H = PR_Hess(z,A,b);
    Hess_Time = toc(Hess_Start);
    real_H = real(4*J'*H*J);
    
    Eig_Prob_Start = tic;
    
    % Call the subroutine CR_Sub_QEP to solve the subproblem
    [d, t, sub_k, sub_tol, GEP_Sol, prep_time, PM_time] = CR_Sub_QEP(real_g,real_H,sigma, GEP_Init_Pt,k,sub_tol_mode);
    
    Eig_Prob_Time = toc(Eig_Prob_Start);
    
    % Warm-Start of GEP: Using current iteration's solution as the initial point
    % for the next iteration's subproblem
    GEP_Init_Pt = GEP_Sol;
    
    x_new = [real(z); imag(z)] + d;
    z = x_new(1:n,1) + 1i*x_new(n+1:2*n,1);
    val = PR_f(z,A,b);
    rel_err = PR_dist(z,z_true)/norm(z_true);
    err_vec(k) = rel_err;
    
%     % system err refers to the error of the optimality system:
%     %      (H + sigma t I)d = -g,   t = |d|,   H + sigma t I >= 0 .
%     sys_err_1 = norm((real_H+sigma*eye(2*n)*t)*d - real_g);
%     sys_err_2 = abs(t - norm(d));
    
    Subprob_Time = toc(Subprob_Start);

    fprintf('k = %d:  Rel. Error = %f,  Func. Val. = %f, Subprob. Tol = %f, Total Subprob. Time = %f\n', k, rel_err, val, sub_tol, Subprob_Time);
    fprintf('Subprob. Details: Grad. = %f,  Hess. = %f,  Eig. Prob. = %f\n', Grad_Time, Hess_Time, Eig_Prob_Time);
    fprintf('SIPM Details: Iter. No. = %d, rightmost eig = %f, Prep. Time = %f,  Loop Time = %f,  Loop Time/Iter. = %f\n', sub_k, t, prep_time, PM_time, PM_time/sub_k);
%     fprintf('system_err_1 = %f, system_err_2 = %f\n\n', sys_err_1, sys_err_2);
end


overall_time = toc(start_time);
fprintf('Init. Time = %f,  Total Time = %f.\n',Init_Time,overall_time);
figure, semilogy(err_vec); 
ylabel('log(Relative Error)');
xlabel('Iteration');