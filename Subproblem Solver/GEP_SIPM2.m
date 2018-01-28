function [x, lambda, k, prep_time, PM_time] = GEP_SIPM2(g,H,sigma,tol,k_max,x_Init,s)

% This is an implementation of the SIPM for solving the inner problem of
% cubic regularization applied to phase retrieval.

n = length(g);
n_hat = 2*n+1;

prep_start = tic;

% Below are the linearization step and the preparation of some quantities that will be used in
% the power iterations.
N1 = sparse(blkdiag(1/s, (1/sigma)*fliplr(eye(2*n)))); % Inverse of (C + s*M)
R = [0, zeros(1,n), (1/s)*g'; (1/sigma)*g, (1/sigma)*H, zeros(n,n); zeros(n,1), (-1/sigma)*eye(n), (1/sigma)*H];
s1_s = tic;
[L2,U2,P2] = lu(s*eye(n_hat) + R); s1 = toc(s1_s); s2_s = tic;
N2 = U2\(L2\P2); s2 = toc(s2_s); s3_s = tic;
n1 = N1(:,1);
q1 = N2*n1; %N2*N1(:,1);
s3 = toc(s3_s);
s4_s = tic;
Q2 = R*N2; s4 = toc(s4_s); s5_s = tic;
q3 = R*q1;
s5 = toc(s5_s);
% fprintf('s1 = %f, s2 = %f, s3 = %f, s4 = %f, s5 = %f\n',s1,s2,s3,s4,s5);
prep_time = toc(prep_start);

k = 0;
x = x_Init;
acc = Inf;
PM_start = tic;
while ((acc > tol) && (k <= k_max))
    % fprintf('sub_k = %d,  sub_acc = %8.2E\n',k,acc);
    x1 = x(1:n_hat,1);
    x2 = x(n_hat+1:2*n_hat,1);
    k = k + 1;
    y1 = N2*x1 + x2(1)*q1;
    y2 = x2(1)*n1 - Q2*x1 - x2(1)*q3; %N1*x2 - Q2*x1 - x2(1)*q3;
    y = [y1; y2];
    norm_y = norm(y);
    mu = x'*y;
    acc = norm(x - y./mu)/norm_y;
    x = y/norm_y;
end
PM_time = toc(PM_start);
% fprintf('Prep. time = %f, PM time = %f, PM time/Iter = %f\n', prep_time, PM_time, PM_time/k);
lambda = s + (1/mu);