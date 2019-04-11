addpath('C:\Users\mayn\Documents\MATLAB\NL0_sphere\solver')
addpath('C:\Users\mayn\Documents\MATLAB\NL0_sphere\data')
addpath('C:\Users\mayn\Documents\MATLAB\lowrankquadrmax-master\matlab\bin')

%% *********** read the file *****************************

file_path = 'C:\Users\mayn\Documents\MATLAB\NL0_sphere\data\face\';

all_file_data = zeros(2429, 361);
pgm_path_list = dir(strcat(file_path, '*.pgm')); % read all file in file_path 
img_num = length(pgm_path_list); % calculate the number of file

for i = 1:img_num
    file_name = strcat(file_path, pgm_path_list(i).name); % generate the name of each file
    file_data = reshape(imread(file_name), 1, []); % read the file and reshape it
    all_file_data(i,:) = file_data; 
end

%% ********** PGM ****************
samplemean = mean(all_file_data, 1);
data = all_file_data - repmat(samplemean, 2429, 1);
A = data' * data /2428;

n=size(A,1);
options.tol = 1e-6;
options.issym = 1;
options.disp  = 0;
options.v0 = ones(n,1);
[xint,Asnorm] =eigs(@(y)(A*y),n,1,'LM',options);

OPTIONS_PGM.tol = 1.0e-6;
OPTIONS_PGM.printyes = 0;
OPTIONS_PGM.Lipconst = 3*Asnorm;
OPTIONS_PGM.maxiter = 3000;
gamma = 1.0e-5;

[xopt,iter,variance] = PGM_PosL0_sphere(abs(xint), -A, OPTIONS_PGM,gamma, 40);

% xopt = xopt - repmat(mean(xopt,1), n, 1);

new_x = reshape(xopt, 19, 19);
figure(1)
imshow(mat2gray(new_x))

%% *********** spannc ***********************

params.algorithm     = 'nnsparse';        % Seek sparse solution
params.nnz           = 40;                % Target sparsity values
params.apprxrank     = 3;                 % Rank of approximation
params.maxsamples    = 2429;               % Number of subspace samples

% (Optional)
params.inputdata     = 'rows';            % Samples in Y are rows
params.maxnoupditer  = 1e3;               % Max iters without improvement    
params.centerdata    = true;              % Center samples (subtract mean)
params.standardata   = false;             % Standardize features
params.logfile       = 'sample.log';      % Path to log file
params.logfilelevel  = 'off';             % Log level for log file
params.logcwlevel    = 'info';            % Log level for command window

[X] = spanpc(all_file_data, params);
xopt = full(X);
%xopt = xopt - repmat(mean(xopt,1), n, 1);
new_x = reshape(xopt, 19, 19);
figure(2)
imshow(mat2gray(new_x))


    
    
    