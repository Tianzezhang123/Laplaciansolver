function [mu] = solveieqn(kernel,alpha,f)
% solve the integral equation [alpha*mu(x) + \int kernel(x,y)mu(y)dy =
% f(x)]
D = diag(alpha);
A = D + kernel;
mu = A\f;
end

