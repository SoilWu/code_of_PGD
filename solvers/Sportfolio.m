%% *************************************************************
% Filename: Sportfolio
%% ************************************************************* 
%%
%% ew means estimate window, kappa means the sparsity of portfolio
%%
%% ***************************************************************
%% Copyright by Yuqia Wu and Shaohua Pan: 2019/4/2
%% ***************************************************************
%% 

function [var, SR, turnover, x_bottle] = Sportfolio(file_name, ew, kappa)

ret = xlsread(file_name);

[nr, nc] = size(ret);

n = fix(nr/5);   % generate weekly data

wk_ret = ret(1:5:1+5*(n-1),:);

ret_rate = (wk_ret(2:n,:)-wk_ret(1:n-1,:))./wk_ret(1:n-1,:); 

n = n-1;        % Number of rows in ret_rate

x_bottle = [];

for i = 1:n-ew
    
    Xw = ret_rate(i:i+ew-1,:);   % estimation window
    
    A_t = Xw'*Xw ; % estimation window covariance matrix
    
    A_t = A_t/norm(A_t,'fro');
    
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

var = sum((sum(ret_rate(ew+1:n,:).*x_bottle(:,1:n-ew)',2)-mu).^2)/(n-ew-1);

sta_d = sqrt(var);    % Standard deviation

SR =  mu/sta_d;       % Sharpe ratio

turnover = sum(sum(abs(x_bottle(:,1:n-ew-1)-x_bottle(:,2:n-ew))));

turnover = turnover/(n-ew-1);