function [A,b,z_true] = PR_Inst(m,n,noise_mag)
% generate the PR problem

if (nargin < 3)
    noise_mag = 0;
end

z_true = (1/sqrt(2))*(randn(n,1) + 1i*randn(n,1));

A = (1/sqrt(2))*(randn(n,m) + 1i*randn(n,m));

b = abs(A'*z_true) + noise_mag*randn(m,1);