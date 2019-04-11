%% *************************************************************
% Filename: APG_L0simplex
%% ************************************************************* 
% Solve the following nonconvex nonsmooth problem with accerlated projection gradient method
% 
%    min f(x)  s.t. ||x||_0<=kappa,||x||_1=1, x>=0                       £¨*£©
%
% In each step of PGM, the following optimization problem is solved:
%
%   min 0.5||x-gk||^2  s.t. ||x||_0<=kappa, ||x||_1=1, x>=0              (**)
% 
% with gk = yk - tau*grad f(yk) and the problem (**) is solved by PosL0_sphere
%% *************************************************************
%%
%% ***************************************************************
%% Copyright by Yuqia Wu and Shaohua Pan: 2018/03/10
%% ***************************************************************
%% 


function [xopt,iter,objval] = APGD_L0simplex(x,A,OPTIONS,gamma,kappa)
%%
if isfield(OPTIONS,'maxiter');    maxiter    = OPTIONS.maxiter;    end
if isfield(OPTIONS,'printyes');   printyes   = OPTIONS.printyes;   end
if isfield(OPTIONS,'Lipconst');   Lipconst   = OPTIONS.Lipconst;   end
if isfield(OPTIONS,'tol');        tol        = OPTIONS.tol;        end

%%
if (printyes)
    fprintf('\n *****************************************************');
    fprintf('******************************************');
    fprintf('\n ************* PGM for the nonnegative sparse sphere constraint **************');
    fprintf('\n ****************************************************');
    fprintf('*******************************************');
    fprintf('\n  iter   optmeasure     diff_obj     norm_gdiff     gamma      tau    time');
end

%% ***************** set the step-size tau ***********************

gamma_Lip = gamma+Lipconst;

tau = 1/gamma_Lip;

%%
%% ***********************  Main Loop ***************************

tstart = clock;

xold = x;

Ax = A*x;  Axold = Ax; 

loss_old = 0;

n = size(A, 1);

for iter = 1:maxiter
    
    if (iter == 1)  % Notice: steepest descent step is the key!!!
        
        t = 1;
        
        yk = x;
        
        gradyk = 2*Ax;
       
    else
        
        t = 0.5*(1+sqrt(1+4*told^2));
        
        beta = (told-1)/t;
        
        beta1 = (1+beta);
        
        beta2 = beta;
        
        yk = beta1*x - beta2*xold;
        
        Ayk = beta1*Ax - beta2*Axold;
         
        gradyk = 2*Ayk;
     
    end
    
    gk = yk - tau*gradyk;
    
    [sg,idx] = sort(gk,'descend'); 

    xk_proj = Proj_simplex(sg(1:kappa),1);
    
    xnew = zeros(n, 1);
    
    xnew(idx(1:kappa)) = xk_proj;
    
    Axnew = A*xnew;
    
    loss = xnew'*Axnew;
    
    ttime = etime(clock,tstart);
    
   %% ************** generate the new iterate xnew ****************
    
    xdiff = xnew - yk;
         
    if (iter == 1)
         grad_diff = Axnew - Ax;
    else    
        
         grad_diff = Axnew-Ayk;
    end
    
    residual = 2*grad_diff - 2*gamma_Lip*xdiff;
    
    normx = norm(xnew);
    
    opt_measure = norm(residual)/max(1,normx);
    
    norm_gdiff = norm(grad_diff)/max(1,normx);
    
    diff_obj = loss - loss_old;
    
    if (printyes)&&(mod(iter,1)==0)
        
        fprintf('\n %3d    %3.2e     %3.2e     %3.2e      %3.2e     %3.2e    %3.1f',iter,opt_measure,diff_obj,norm_gdiff,gamma,tau,ttime);
        
    end
    %%
    %% ************* check stopping criterion ******************
    %%
    if (opt_measure<tol)
       
        xopt = xnew;
        
        objval = loss;
        
        return;
    end
    
    xold = x;  x = xnew; 
    
    Axold = Ax; Ax = Axnew;  
    
    loss_old = loss;
    
     told = t;
end

if (iter==maxiter)
    
    xopt = x;
    
   objval = loss;
    
end


