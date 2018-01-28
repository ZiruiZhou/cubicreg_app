function err = LRM_dist(U,U_true)
% distance of U to the optimal solution set

U_tTU = U_true'*U;
[P,S,R] = svd(U_tTU);
err = norm(U - U_true*P*R');
