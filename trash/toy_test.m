%% ********************** generate the data ***********************************

v1 = [1, 1, 1, 1, 0, 0, 0, 0, 0.9, 0.9]';
v2 = [0, 0, 0, 0, 1, 1, 1, 1, -0.3, 0.3]';

v1 = v1/ norm(v1, 2);
v2 = v2/ norm(v2, 2);

n = size(v1, 1);
I = eye(n);
P = randn(n, n);
P(:, 1:2) = [v1, v2];
for i = 3:10
    P(:, i) = P(:, i)/norm(P(:, i));
end

P = Gram_Schmidt(P);

d = [200, 100, 50, 50, 6, 5, 4, 3, 2, 1];
D = diag(d);

Sigma = P *D *P';
mu = zeros(1,n);
X = mvnrnd(mu,Sigma,100);
A = X'*X;

%% *************************** NPGM *************************************

n=size(A,1);
options.tol = 1e-6;
options.issym = 1;
options.disp  = 0;
options.v0 = ones(n,1);
[xint,Asnorm] =eigs(@(y)(A*y),n,1,'LM',options);

OPTIONS_PGM.tol = 1.0e-6;
OPTIONS_PGM.printyes = 0;
OPTIONS_PGM.Lipconst = Asnorm;
OPTIONS_PGM.maxiter = 3000;
gamma = 1.0e-5;

[xopt1,iter, ~] = PGM_PosL0_sphere(abs(xint), -A, OPTIONS_PGM,gamma, 4);



A1 = (I - xopt1 * xopt1') * A * (I - xopt1 * xopt1');
n=size(A,1);
options.tol = 1e-6;
options.issym = 1;
options.disp  = 0;
options.v0 = ones(n,1);
[xint,Asnorm] =eigs(@(y)(A1*y),n,1,'LM',options);

OPTIONS_PGM.tol = 1.0e-6;
OPTIONS_PGM.printyes = 0;
OPTIONS_PGM.Lipconst = Asnorm;
OPTIONS_PGM.maxiter = 3000;
gamma = 1.0e-5;

[xopt2,iter, ~] = PGM_PosL0_sphere(abs(xint), -A1, OPTIONS_PGM,gamma, 4);

%% ********************** Spnn ******************************

params.algorithm     = 'nnsparse';        % Seek sparse solution
params.nnz           = 4;                % Target sparsity values
params.apprxrank     = 3;                 % Rank of approximation
params.maxsamples    = 100;               % Number of subspace samples

% (Optional)
params.inputdata     = 'rows';            % Samples in Y are rows
params.maxnoupditer  = 1e3;               % Max iters without improvement    
params.centerdata    = true;              % Center samples (subtract mean)
params.standardata   = false;             % Standardize features
params.logfile       = 'sample.log';      % Path to log file
params.logfilelevel  = 'off';             % Log level for log file
params.logcwlevel    = 'info';            % Log level for command window

[K] = spanpc(X, params);
xopt3 = full(K);

X = (I - xopt3 * xopt3') * X';

params.algorithm     = 'nnsparse';        % Seek sparse solution
params.nnz           = 4;                % Target sparsity values
params.apprxrank     = 3;                 % Rank of approximation
params.maxsamples    = 100;               % Number of subspace samples

% (Optional)
params.inputdata     = 'columns';            % Samples in Y are rows
params.maxnoupditer  = 1e3;               % Max iters without improvement    
params.centerdata    = true;              % Center samples (subtract mean)
params.standardata   = false;             % Standardize features
params.logfile       = 'sample.log';      % Path to log file
params.logfilelevel  = 'off';             % Log level for log file
params.logcwlevel    = 'info';            % Log level for command window

[K] = spanpc(X, params);
xopt4 = full(K);


xopt1' * A * xopt1
xopt2' * A * xopt2

xopt3' * A * xopt3
xopt4' * A * xopt4