function [xopt] = Spnn_solver(X, k) 


params.algorithm     = 'nnsparse';        % Seek sparse solution
params.nnz           = k;                % Target sparsity values
params.apprxrank     = 3;                 % Rank of approximation
params.maxsamples    = 165;               % Number of subspace samples

% (Optional)
params.inputdata     = 'columns';            % Samples in Y are rows
params.maxnoupditer  = 1e3;               % Max iters without improvement    
params.centerdata    = true;              % Center samples (subtract mean)
params.standardata   = false;             % Standardize features
params.logfile       = 'sample.log';      % Path to log file
params.logfilelevel  = 'off';             % Log level for log file
params.logcwlevel    = 'info';            % Log level for command window

[Y] = spanpc(X, params);
xopt = full(Y);
%xopt = xopt - repmat(mean(xopt,1), n, 1);
% new_x = reshape(xopt, 100, 100);