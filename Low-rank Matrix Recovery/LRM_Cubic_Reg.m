function [U, rel_err, val, overall_time] = LRM_Cubic_Reg(U_true, A, b, init, tol, k_max)

% This program solves the low-rank matrix recovery problem by cubic regularization.
% The inner problem is handled as a QEP (quadratic eigenvalue problem).

if (nargin == 3)
    init = 0;  % default is random arbitrary initialization
    tol = 10^(-4);
    k_max = 100;
elseif (nargin == 4)
    tol = 10^(-4);
    k_max = 100;
elseif (nargin == 5)
    k_max = 1000;
end

[n,r] = size(U_true);
m = size(A,3);
nr = n*r;

% Initialization
Init_Start = tic;
if init==0
    U0 = 5*randn(n,r);
else
    fprintf('Initializing...,');
    U0 = LRM_Init(A,b,m,n,r);
    fprintf('OK! Finished. Time consumed: %f\n\n\n', toc(Init_Start));
end
Init_Time = toc(Init_Start);

% Lipschitz constant
L2 = 10*norm(U0);

% Initial Point of SIPM for subproblem
GEP_Init_Pt = randn(4*nr+2,1);
GEP_Init_Pt = GEP_Init_Pt/norm(GEP_Init_Pt);

% Subproblem error tolerance mode: 
%     1 stands for constant tolerance; 
%     2 stands for diminising tolerance.
sub_tol_mode = 2;

start_time = tic;
k = 0;
rel_err = Inf;
U = U0;

% Evaluate the function value, gradient, and Hessian.
Evaluation_Start = tic;
[val,g,H] = LRM_f(U,A,b);
Evaluation_Time = toc(Evaluation_Start);

err_vec = zeros(k_max,1);
while ((rel_err > tol) && (k <= k_max))
    Subprob_Start = tic;
    k = k + 1;
    sigma = L2;
    
    
    
        
    % Call the subroutine CR_Sub_QEP to solve the subproblem.
    Eig_Prob_Start = tic;
    [d, t, sub_k, sub_tol, GEP_Sol, prep_time, PM_time] = CR_Sub_QEP(g,H,sigma, GEP_Init_Pt,k,sub_tol_mode);
    Eig_Prob_Time = toc(Eig_Prob_Start);
    
    
    % Warm-Start of GEP: Using current iteration's solution as the initial point
    % for the next iteration's subproblem
    GEP_Init_Pt = GEP_Sol;
    
    U = U + reshape(d,n,r);
    
    Evaluation_Start = tic;
    [val,g,H] = LRM_f(U,A,b);
    Evaluation_Time = toc(Evaluation_Start);
    
    rel_err = LRM_dist(U,U_true)/norm(U_true);
    err_vec(k) = rel_err;
    
%     % system err refers to the error of the optimality system:
%     %      (H + sigma t I)d = -g,   t = |d|,   H + sigma t I >= 0 .
%     sys_err_1 = norm((real_H+sigma*eye(2*n)*t)*d - real_g);
%     sys_err_2 = abs(t - norm(d));
    
    Subprob_Time = toc(Subprob_Start);

    fprintf('k = %d:  Rel. Error = %f,  Func. Val. = %f, Subprob. Tol = %f, Total Subprob. Time = %f\n', k, rel_err, val, sub_tol, Subprob_Time);
    fprintf('Subprob. Details: Evaluation Time = %f,  Eig_Prob_Time = %f\n', Evaluation_Time, Eig_Prob_Time);
    fprintf('SIPM Details: Iter. No. = %d, rightmost eig = %f, Prep. Time = %f,  Loop Time = %f,  Loop Time/Iter. = %f\n', sub_k, t, prep_time, PM_time, PM_time/sub_k);
%     fprintf('system_err_1 = %f, system_err_2 = %f\n\n', sys_err_1, sys_err_2);
end

overall_time = toc(start_time);
fprintf('Init. Time = %f,  Total Time = %f.\n',Init_Time,overall_time);
figure, semilogy(err_vec); 
ylabel('log(Relative Error)');
xlabel('Iteration');