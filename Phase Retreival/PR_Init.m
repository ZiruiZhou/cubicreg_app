function z0 = PR_Init(A,b)
% Initialization for PR problem

n = size(A,1);

npower_iter = 100;                           % Number of power iterations 
z0 = randn(n,1); z0 = z0/norm(z0,'fro');    % Initial guess 
for tt = 1:npower_iter,                     % Power iterations
    z0 = A*((b.^2).* (A'*z0));
    z0 = z0/norm(z0,'fro');
end

normest = sqrt(sum(b.^2)/numel(b));    % Estimate norm to scale eigenvector  
z0 = normest * z0;                   % Apply scaling