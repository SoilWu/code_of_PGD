addpath(genpath(pwd));
addpath('C:\Users\mayn\Documents\MATLAB\NL0_sphere\solver')
addpath('C:\Users\mayn\Documents\MATLAB\NL0_sphere\data')
addpath('C:\Users\mayn\Documents\MATLAB\lowrankquadrmax-master\matlab\bin')

%% *********** read the file *****************************

file_name = 'portfolio.xlsx';

X_t = xlsread(file_name);

[n, p] = size(X_t);

idx = [];

for i = 1:p
    if X_t(1,i)<-40
        idx = [idx, i];
    end
end

nn = size(idx, 2);

for i=1:nn
    X_t(:,idx(i)-i+1)=[];
end

[n, p] = size(X_t);

X_t = X_t(:, 2:p);

X = zeros(n, p-1);

X(1,:) = X_t(1,:);

for i = 2:n
    
    X(i,:) = X(i-1,:)+X_t(i,:);  % Compute the price
    
end


ew = 120;

kappa = 30;

x_bottle = [];

tic;

for i = 1:n-ew
    
    X_w = X(i:i+119,:); % estimation window
    
    A_t = X_w' * X_w ; % estimation window covariance matrix

    A_t = A_t/norm(A_t, 2);
    
    [xopt, iter] = Sportfolio (A_t, kappa);

    x_bottle(:,i) = xopt;
    
    fprintf('\n %3d  %3d',i,iter);
    
end

toc

mu = trace(x_bottle'*X_t(ew+1:n,:)')/(n-ew);
sigma = 0;
for i = 1:991
    kk = X_t(ew+i,:)* x_bottle(:,i);
    sigma = sigma + (kk - mu)^2;
end 
sigma = sigma / 990;
SR =  mu/sqrt(sigma);
turnover = 0;

for i = 1:990
    turnover = turnover + sum(abs(x_bottle(:,i)-x_bottle(:,i+1)));
end
turnover = turnover/991;
    
    
    

