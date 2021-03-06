%% *************************************************************
% Filename: PGD_PosL0_sphere
%% ************************************************************* 
% Solve the following nonconvex nonsmooth problem with projection gradient method
% 
%    min f(x)  s.t. ||x||_0<=kappa,||x||=1, x>=0                       ��*��
%
% In each step of PGM, the following optimization problem is solved:
%
%   min 0.5||x-gk||^2  s.t. ||x||_0<=kappa, ||x||=1, x>=0              (**)
% 
% with gk=x^k - tau*grad f(x^k), and the problem (**) is solved by PosL0_sphere
%% *************************************************************
%%
%% ***************************************************************
%% Copyright by Yuqia Wu and Shaohua Pan: 2018/12/23
%% ***************************************************************
%% 

function [xopt,iter,svariance] = PGD_PosL0_sphere(x,A,OPTIONS,gamma,kappa)
%%
if isfield(OPTIONS,'maxiter');    maxiter    = OPTIONS.maxiter;    end
if isfield(OPTIONS,'printyes');   printyes   = OPTIONS.printyes;   end
if isfield(OPTIONS,'Lipconst');   Lipconst   = OPTIONS.Lipconst;   end
if isfield(OPTIONS,'tol');        tol        = OPTIONS.tol;        end

%%
if (printyes)
    fprintf('\n *****************************************************');
    fprintf('******************************************');
    fprintf('\n ************* PGD for the nonnegative sparse sphere constraint **************');
    fprintf('\n ****************************************************');
    fprintf('*******************************************');
    fprintf('\n  iter   optmeasure     diff_obj     norm_gdiff     gamma      tau    time');
end

%% ***************** set the step-size tau ***********************

gamma_Lip = gamma + Lipconst;

tau = 0.5/gamma_Lip;

%%
%% ***********************  Main Loop ***************************

tstart = clock;

Ax = A*x;  loss_old = 0;

optmeasure_list = zeros(maxiter,1);

for iter = 1:maxiter
    
    gradf = 2*Ax;
    
    gk = x - tau*gradf;
    
    xnew = PosL0_sphere(gk,kappa);
    
    Axnew = A*xnew;
    
    loss = xnew'*Axnew;
    
    ttime = etime(clock,tstart);
    
   %% ************** generate the new iterate xnew ****************
    
    xdiff = xnew - x;
    
    grad_diff = Axnew - Ax;
    
    residual = 2*(grad_diff - gamma_Lip*xdiff);
    
    normx = norm(xnew);
   
    opt_measure = norm(residual)/max(1,normx);
    
    optmeasure_list(iter) = opt_measure;
    
    norm_gdiff = norm(grad_diff)/max(1,normx);
           
    diff_obj = loss - loss_old;
    
    if (printyes)&&(mod(iter,10)==0)
        
        fprintf('\n %3d    %3.2e     %3.2e     %3.2e      %3.2e     %3.2e    %3.1f',iter,opt_measure,diff_obj,norm_gdiff,gamma,tau,ttime);
        
    end
    %%
    %% *************** check stopping criterion ******************
    %%
    if (opt_measure<=tol)||(iter>50 && opt_measure<=3*tol && max(abs(optmeasure_list(iter-50:iter)))<=10*tol) 
        
        xopt = xnew;
        
        svariance = -loss;
        
        return;
   
    end
    
    x = xnew;  Ax = Axnew;  
    
    loss_old = loss;
    
end

if (iter==maxiter)
    
    xopt = xnew;
    
    svariance = -loss;
    
end


