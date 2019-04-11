%% ************ Cumulative variance of the first few PCs *************

addpath(genpath(pwd))

%% *********** Compute the covariance matrix **********************

X = xlsread('\data\leukemia.xlsx');

[p,n] = size(X);

mean_x = mean(X,2);

tempX = X - (mean(X,2)*ones(1,n));

A = tempX*tempX'/71;

TrA = trace(A);

k = 10;             % the number of PCs we will solve

%%
%% **************** sparsity =1000 *******************************

sparsity = 1000;

y1 = zeros(k,1);

y3 = zeros(k,1);

PGD_time = zeros(k,1);

Spnn_time = zeros(k,1);

%********************** Proximal gradient descent method ********************

OPTIONS_PGD.tol = 1.0e-6;

OPTIONS_PGD.printyes = 0;

OPTIONS_PGD.maxiter = 3000;

xopt = zeros(p,1);

Ak1 = A;

Akxopt1 = zeros(p,1);

temp_const1 = xopt'*Akxopt1;

svariance = 0 ;

for i = 1:k
    i
    tstart = clock;
    
    tempAk1 = Akxopt1*xopt';
    
    Ak1 = Ak1 - (tempAk1+tempAk1') + (temp_const1*xopt)*xopt';
    
    options.tol = 1e-8;
    
    options.issym = 1;
    
    options.disp  = 0;
    
    options.v0 = ones(p,1);
    
    [xint,variance] = eigs(@(y)(Ak1*y),p,1,'LM',options);
        
    OPTIONS_PGD.Lipconst = variance;
    
    [xopt,iter] = PGD_PosL0_sphere(abs(xint),-Ak1,OPTIONS_PGD,1e-5,sparsity);
           
    svariance = svariance + xopt'*(A*xopt);
    
    y1(i) = svariance/TrA;
    
    Akxopt1 = Ak1*xopt;
    
    temp_const1 = xopt'*Akxopt1;
    
    PGD_time(i) = etime(clock,tstart);
end

%************************* Spannogram ************************************

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

svariance = 0 ;

X1 = X;

for i = 1:k
    
    tstart = clock;
    
    [K] = spanpc(X1,params);
    
    zopt = full(K);
    
    X1 = X1 - zopt*(zopt'* X1);
    
    svariance = svariance + zopt'*(A*zopt);
    
    y3(i) = svariance /TrA;
    
    Spnn_time(i) = etime(clock,tstart);
end

%% ********************* plot figures *****************************

x1 = [1   2   3   4   5   6   7   8   9  10];

x3 = x1;

subplot(1,2,1);
h=plot(x1,y1,'r-*',x3,y3,'b-o');
set(h,'LineWidth',2) 
xlabel('Number of PCs (k)');   
ylabel('Cumulative explained variance(%)');
set(get(gca,'XLabel'),'FontSize',14);
set(get(gca,'YLabel'),'FontSize',14);
legend('Algorithm 1','Spannogram');
title('\kappa = 1000')
grid on
hold on;

%% **************** sparsity =2000 *******************************

sparsity = 2000;

y2 = zeros(k,1);

y4 = zeros(k,1);

%********************** Proximal gradient descent method ********************

xopt = zeros(p,1);

Ak2 = A;

Akxopt2 = zeros(p,1);

temp_const2 = xopt'*Akxopt2;

svariance = 0 ;

for i = 1:k
   i
    tempAk2 = Akxopt2*xopt';
    
    Ak2 = Ak2 - (tempAk2+tempAk2') + (temp_const2*xopt)*xopt';
    
    options.tol = 1e-8;
    
    options.issym = 1;
    
    options.disp  = 0;
    
    options.v0 = ones(p,1);
    
    [xint,variance] = eigs(@(y)(Ak2*y),p,1,'LM',options);
        
    OPTIONS_PGD.Lipconst = variance;
    
    [xopt,iter] = PGD_PosL0_sphere(abs(xint),-Ak2,OPTIONS_PGD,1e-5,sparsity);
   
    svariance = svariance + xopt'*(A*xopt);
    
    y2(i) = svariance/TrA;
    
    Akxopt2 = Ak2*xopt;
    
    temp_const2 = xopt'*Akxopt2;
    
end

%************************* Spannogram ************************************

params.nnz = sparsity;     % Target sparsity values

svariance = 0 ;

X2 = X;

for i = 1:k
    
    [K] = spanpc(X2,params);
    
    zopt = full(K);
    
    X2 = X2 - zopt*(zopt'*X2);
    
    svariance = svariance + zopt'*(A*zopt);
    
    y4(i) = svariance/TrA;
    
end

%% ************************ plot figures *************************

x2 = [1  2   3   4  5  6  7  8  9  10];

x4 = x2;

subplot(1,2,2);
h=plot(x2,y2,'r-*',x4,y4,'b-o');
set(h,'LineWidth',2) 
xlabel('Number of PCs (k)');   
ylabel('Cumulative explained variance(%)');
set(get(gca,'XLabel'),'FontSize',14);
set(get(gca,'YLabel'),'FontSize',14);
legend('Algorithm 1','Spannogram');
title('\kappa = 2000')
grid on
hold on;

plot([1:10],PGD_time,'r-*',[1:10],Spnn_time,'b-o','LineWidth',2);
