function [d, t, sub_k, sub_tol, GEP_Sol, prep_time, PM_time] = CR_Sub_QEP(g,H,sigma,GEP_Init_Pt,k,mode)

% This function solves the QEP by linearization, a process that turns QEP
% into a GEP (generalized eigenvalue problem).

n = length(g);
s = n;

%Setting the mode of subproblem error tolerance
if (mode == 1)
    sub_tol = 10^(-3);  % constant error tolerance
elseif (mode == 2)
    sub_tol = 1/k^2*10^(-4);  % diminishing error tolerance
end

sub_k_max = 10000;
[GEP_Sol, t, sub_k, prep_time, PM_time] =  GEP_SIPM2(g,H,sigma,sub_tol,sub_k_max,GEP_Init_Pt,s);
d = GEP_Sol(2:n+1)./GEP_Sol(1);