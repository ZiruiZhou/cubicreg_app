function f = PR_f(z,A,b)
m = length(b);
f = (norm(abs(A'*z).^2 - b.^2)^2)/(2*m);