function [A,b,U_true] = LRM_Inst(m,n,r,noise_mag)

% This function generate the true solution U, measurement matrix A, and
% measurements b;
% Each A_i, i=1,...,m, is reshaped to a column vector, and stored as a
% column in A. So A is of dimension n^2*m.

if (nargin < 3)
    noise_mag = 0;
end

% Construct the true solution X_true to be recovered.
U_true = -5 + 10*rand(n,r);
X_true = U_true*U_true';
X_vec = reshape(X_true,n*n,1);

% Construct the measurement symmetric matrices A_1,A_2,...,A_m.
% The tensor A_s stores all the matrices A_i's.
A_ns = randn(n,n,m);
A = 0.5*(A_ns + permute(A_ns,[2,1,3]));

% Calculate the measurement vector b.
A_vec = reshape(A,n*n,m);
b = A_vec'*X_vec;