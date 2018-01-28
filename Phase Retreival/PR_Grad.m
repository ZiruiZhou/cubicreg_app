function g = PR_Grad(z,A,b)
% [n,m] = size(A);
% g = zeros(2*n,1);
% for j=1:m
%     w = (abs(A(:,j)'*z)^2 - b(j)^2)*(A(:,j)'*z)*A(:,j);
%     g = g + [w; conj(w)];
% end
% g = (sqrt(2)/m)*g;

% use vecterized operations
m = size(A,2);
g_upper = 1/m*bsxfun(@times, A, (abs(A'*z).^2 - b.^2)')*A'*z;
g = [g_upper; conj(g_upper)];
