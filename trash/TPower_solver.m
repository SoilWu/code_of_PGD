function [xopt, new_x, iter] = TPower_solver(A, sparsity)

n=size(A,1);
options.tol = 1e-6;
options.issym = 1;
options.disp  = 0;
options.v0 = ones(n,1);
[xint,Asnorm] =eigs(@(y)(A*y),n,1,'LM',options);

OPTIONS_TPower.tol = 1.0e-6;
OPTIONS_TPower.printyes = 0;
OPTIONS_TPower.Lipconst = Asnorm;
OPTIONS_TPower.maxiter = 3000;
gamma = 1.0e-5;

[xopt,iter, ~] = TPower(-abs(xint), -A, OPTIONS_TPower, sparsity);
new_x = reshape(xopt, sqrt(n), sqrt(n));
