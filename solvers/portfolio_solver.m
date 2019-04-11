%% *************************************************************
% Filename: portfolio_solver
%% ************************************************************* 
%%
%% ew means estimate window, kappa means the sparsity of portfolio
%%
%% *************************************************************
%%
%% ***************************************************************
%% Copyright by Yuqia Wu and Shaohua Pan: 2019/4/2
%% ***************************************************************
%% 

function [var, SR, turnover, x_bottle] = portfolio_solver(file_name, ew, kappa)

ret = xlsread(file_name);

[n, nc] = size(ret);

k = fix(n/5); % generate weekly data

wk_ret = [];

for i = 1:k
    
    wk_ret(i,:) = ret(1+5*(i-1),:);
    
end

[n, p] = size(wk_ret);

ret_rate = [];
    
for i = 2:n
    
    ret_rate(i-1, :) = (wk_ret(i, :) - wk_ret(i-1, :))./wk_ret(i-1,:); % rate of return
    
end

n = n-1; %  Number of rows in ret_rate

x_bottle = [];

for i = 1:n-ew
    
    X_w = ret_rate(i:i+ew-1,:); % estimation window
    
    A_t = X_w' * X_w ; % estimation window covariance matrix

    A_t = A_t/norm(A_t, 'fro');
    
     options.tol = 1e-6;
    options.issym = 1;
    options.disp  = 0;
    options.v0 = ones(nc,1);  % nc is the number of columns of A_t
    
    [xint,Asnorm] =eigs(@(y)(A_t*y),nc,1,'LM',options);
    
    OPTIONS_PGD.tol = 1.0e-6;
    
    OPTIONS_PGD.printyes = 0;
    
    OPTIONS_PGD.Lipconst = Asnorm;
    
    OPTIONS_PGD.maxiter = 30000;
    
    gamma = 1.0e-5;
    
    xopt = PGD_L0simplex(abs(xint),A_t,OPTIONS_PGD,gamma,kappa);

    x_bottle(:,i) = xopt;
end

mu = trace(x_bottle'*ret_rate(ew+1:n,:)')/(n-ew);

var = 0;

for i = 1:n-ew
    kk = ret_rate(ew+i,:)* x_bottle(:,i);
    var = var + (kk - mu)^2;
end 
var = var / (n-ew-1); % var
sta_d = sqrt(var); % Standard deviation
SR =  mu/sta_d; % sharp ratio
turnover = 0;

for i = 1:n-ew-1
    turnover = turnover + sum(abs(x_bottle(:,i)-x_bottle(:,i+1))); 
end
turnover = turnover/(n-ew-1);