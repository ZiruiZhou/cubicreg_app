function [rel_err,overall_time] = LRM_Test(n,r,sav)

% LRM_Test(n,sav) test the performance of the cubic regularization method
% applied to Low-rank Matrix Recovery problem.

if (nargin<2)
    error('Error of Input: Users need to specify both n and r.');
elseif (nargin<3)
    sav = 0;   % do not save the data of this run
end

% For simple test, use the following one-dimensional problem
%     n = 1;
%     z_true = 1;
%     A = [1,1i];
%     b = abs(A'*z_true);

addpath('..\Subproblem Solver');


m = ceil(3*n*r);   % number of measurements

[A,b,U_true] = LRM_Inst(m,n,r);    % generate the problem with random matrices
% A is a tensor with each slice is a symmetric measurement matrix A_i.

[~, rel_err, ~, overall_time] = LRM_Cubic_Reg(U_true, A, b,0,10^(-8),100);

if sav
    filename = sprintf('%s_%d_%d','CR_LRM',n,r);
    save(filename, 'A','b','U_true','n','r','rel_err','overall_time');
end