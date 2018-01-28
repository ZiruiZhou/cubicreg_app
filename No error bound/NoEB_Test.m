% This script test the convergence rate of CR method applied to solve
%               min f(x) = 1/4*x^4.

addpath('..\Subproblem Solver');
sav = 1;

tol = 10^(-6);
k_max = 250;
n = 1;
x_opt = 0;

L2 = 2;

GEP_Init_Pt = randn(4*n+2, 1);
GEP_Init_Pt = GEP_Init_Pt/norm(GEP_Init_Pt);

sub_tol_mode = 2;

start_time = tic;
k = 0;
rel_err = Inf;
x = 1;


val = x^4/4;
g = x^3;
H = 3*x^2;

err_vec = zeros(k_max,1);

while ((rel_err > tol) && (k <= k_max))
    k = k + 1;
    sigma = L2;
    
    Eig_Prob_Start = tic;
    [d,t,sub_k,sub_tol,GEP_Sol,prep_time,PM_time] = CR_Sub_QEP(g,H,sigma, GEP_Init_Pt,k,sub_tol_mode);
    Eig_Prob_Time = toc(Eig_Prob_Start);
    
    GEP_Init_Pt = GEP_Sol;
    
    x = x + d;
    
    val = x^4/4;
    g = x^3;
    H = 3*x^2;
    
    rel_err = norm(x);
    err_vec(k) = rel_err;
    
    fprintf('k = %d:  Rel. Error = %f,  Func. Val. = %f\n', k, rel_err, val);
    fprintf('Subprob. Details: Eig_Prob_Time = %f\n', Eig_Prob_Time);
    fprintf('SIPM Details: Iter. No. = %d, rightmost eig = %f, Prep. Time = %f,  Loop Time = %f,  Loop Time/Iter. = %f\n', sub_k, t, prep_time, PM_time, PM_time/sub_k);
end

overall_time = toc(start_time);
fprintf('Total Time = %f.\n',overall_time);
figure, semilogy(err_vec); 
ylabel('log(Relative Error)');
xlabel('Iteration');


if sav
    save('no_error_bd', 'rel_err','overall_time');
end