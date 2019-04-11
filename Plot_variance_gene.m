addpath(genpath(pwd))

%% *************** Compute the covariance matrix **********************

X = xlsread('\data\leukemia.xlsx');

[p,n] = size(X);

mean_x = mean(X,2);

tempX = X-(mean(X,2)*ones(1,n));

A = tempX*tempX'/71;

%% ************** sparse covariance under different sparsity ***********

sparsity = linspace(100,8100,41);

nn = size(sparsity,2);

x1 = zeros(nn,1);

y1 = zeros(nn,1);

x2 = zeros(nn,1);

y2 = zeros(nn,1);

%********************** Proximal gradient descent method ******************** 

options.tol = 1e-8;

options.issym = 1;

options.disp = 0;

options.v0 = ones(p,1);

[xint,variance] = eigs(@(y)(A*y),p,1,'LM',options);

OPTIONS_PGD.tol = 1.0e-6;

OPTIONS_PGD.printyes = 0;

OPTIONS_PGD.Lipconst = variance;

OPTIONS_PGD.maxiter = 3000;

for i = 1:nn
    i
    [xopt,iter,svariance] = PGD_PosL0_sphere(abs(xint),-A,OPTIONS_PGD,1e-5,sparsity(i));
    
    x1(i) = sparsity(i);
    
    y1(i) = svariance/variance;
    
end

% ************************* Spannogram ************************************ 

params.algorithm     = 'nnsparse';        % Seek sparse solution
params.nnz           = sparsity ;         % Target sparsity values
params.apprxrank     = 3;                 % Rank of approximation
params.maxsamples    = 72;               % Number of subspace samples
params.inputdata     = 'columns';         % Samples in Y are rows
params.maxnoupditer  = 1e3;               % Max iters without improvement    
params.centerdata    = true;              % Center samples (subtract mean)
params.standardata   = false;             % Standardize features
params.logfile       = 'sample.log';      % Path to log file
params.logfilelevel  = 'off';             % Log level for log file
params.logcwlevel    = 'info';            % Log level for command window

[X] = spanpc(X,params);

X = full(X);

for i = 1:nn
           
    x2(i) = sparsity(i);
    
    y2(i) = X(:,i)'*A*X(:,i)/variance;
    
end

plot(x1,y1,'r-+',x2,y2,'b-','LineWidth',2);

xlabel('Sparsity (\kappa)');

ylabel('Proportion of explained variance(%)');

legend('Algorithm 1','Spannogram');

