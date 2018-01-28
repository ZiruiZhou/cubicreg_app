function U0 = LRM_Init(A,b,m,n,r)

% Initialization for LRM problem

% Calculate A^*(b)
A_v = reshape(A,n*n,m);
A_vb = A_v*b;
Atb = reshape(A_vb,n,n);

[U,E] = eigs(Atb,r,'la');
E(E<0) = 0;
E = E.^(0.5);
U0 = U*E;