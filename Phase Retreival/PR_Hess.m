function H = PR_Hess(z,A,b)
m = length(b);

AA = blkdiag(A,conj(A));
D1 = sparse(diag(2*abs(A'*z).^2 - b.^2));
D2 = sparse(diag((A'*z).^2));
D = sparse([D1, D2; conj(D2), D1]);
H = (1/m)*(AA*D)*AA';

% Method 2
% UL = bsxfun(@times, A, (2*abs(A'*z).^2 - b.^2)')*A';
% UR = bsxfun(@times, A, ((A'*z).^2)')*conj(A');
% LL = UR';
% LR = conj(UL);
% H = (2/m)*[UL,UR;LL,LR];

