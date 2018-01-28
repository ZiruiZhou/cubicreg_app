function [rel_err,overall_time] = PR_Test(n,sav)

% PR_Test(n,sav) test the performance of the cubic regularization method
% applied to Phase Retreival problem.

if (nargin<2)
    sav = 0;   % do not save the data of this run
end

% For simple test, use the following one-dimensional problem
%     n = 1;
%     z_true = 1;
%     A = [1,1i];
%     b = abs(A'*z_true);

addpath('..\Subproblem Solver');


m = ceil(3*n*log(n));   % number of measurements

[A,b,z_true] = PR_Inst(m,n);    % generate the problem with random matrices


[~, rel_err, ~, overall_time] = PR_Cubic_Reg(z_true, A, b,0,10^(-8),100);

if sav
    filename = sprintf('%s_%d','CR_PR',n);
    save(filename, 'A','b','z_true','n','rel_err','overall_time');
end