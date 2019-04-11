addpath(genpath(pwd));
addpath('C:\Users\mayn\Documents\MATLAB\NL0_sphere\solver')
addpath('C:\Users\mayn\Documents\MATLAB\NL0_sphere\data')
addpath('C:\Users\mayn\Documents\MATLAB\lowrankquadrmax-master\matlab\bin')

%% *********** read the file *****************************

file_name1 = 'return.xlsx';
file_name2 = 'price.xlsx';

ret = xlsread(file_name1);
price = xlsread(file_name2);

[n, p] = size(ret);

ew = 52;

kappa = 300;

x_bottle = [];

tic;

for i = 1:n-ew
    
    X_w = ret(i:i+ew-1,:); % estimation window
    
    A_t = X_w' * X_w ; % estimation window covariance matrix

    A_t = A_t/norm(A_t, 'fro');
    
    [xopt, iter] = Sportfolio (A_t, kappa);

    x_bottle(:,i) = xopt;
    
    fprintf('\n %3d  %3d',i,iter);
    
end

toc

mu = trace(x_bottle'*ret(ew+1:n,:)')/(n-ew);
sigma = 0;
for i = 1:n-ew
    kk = ret(ew+i,:)* x_bottle(:,i);
    sigma = sigma + (kk - mu)^2;
end 
sigma = sigma / (n-ew-1);
sta_d = sqrt(sigma);
SR =  mu/sta_d;
turnover = 0;

for i = 1:n-ew-1
    turnover = turnover + sum(abs(x_bottle(:,i)-x_bottle(:,i+1)));
end
turnover = turnover/(n-ew-1);