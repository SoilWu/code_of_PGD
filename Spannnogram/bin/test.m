


params.algorithm     = 'nnsparse';        % Seek sparse solution
params.nnz           = [10, 30,  50, 70]; % Target sparsity values
params.apprxrank     = 3;                 % Rank of approximation
params.maxsamples    = 1e4;               % Number of subspace samples

% (Optional)
params.inputdata     = 'rows';            % Samples in Y are rows
params.maxnoupditer  = 1e3;               % Max iters without improvement    
params.centerdata    = true;              % Center samples (subtract mean)
params.standardata   = false;             % Standardize features
params.logfile       = 'sample.log';      % Path to log file
params.logfilelevel  = 'off';             % Log level for log file
params.logcwlevel    = 'info';            % Log level for command window


% Run
[X] = spanpc(all_file_data, params);

% Plot explained variance for each sparsity value
figure;
plot(params.nnz, var(Y*X),'--sr');
title('Explained variance: nonnegative k-sparse principal component')
xlabel('Sparsity (k)');
ylabel('Explained (empirical) Variance');
grid on;
