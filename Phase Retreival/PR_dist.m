function err = PR_dist(z,z_true)
% err = \| z - z_true exp(i*\phi)\|, where \phi = \argmin_[0,2pi] \| z - z_true exp(i*\phi)\|

phi = -angle(z'*z_true);
err = norm(z - z_true*exp(1i*phi));