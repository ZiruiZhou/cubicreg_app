function [f,g,H] = LRM_f(U,A,b)

% Evaluate the function value, gradient, and Hessian at the given point U

m = length(b);
[n,r] = size(U);
Av = reshape(A,n*n,m);
AvU = Av'*reshape(U*U',n*n,1);
diff = AvU - b;

% function value
f = norm(diff)^2/(4*m);

% gradient
comp = reshape(U'*reshape(A,n,n*m), r,n,m);
comp = reshape(permute(comp,[2,1,3]),n*r,m);
g = comp*diff/m;

% hessian
Avd = Av*diff;
AvdCell = repmat({reshape(Avd,n,n)}, 1, r);
H1 = blkdiag(AvdCell{:});
H2 = comp*comp';
H = (H1 + 2*H2)/m;


% Old approaches: (kron product cost more time)
% Ir = eye(r);
% In = eye(n);
% g = kron(U',In)*Avd/m;
% H = (kron(Ir,reshape(Avd,n,n)) + 2*kron(U',In)*Av*Av'*kron(U',In)')/m;



