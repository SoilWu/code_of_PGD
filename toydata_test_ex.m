clear;

%% *********** read the file *****************************

addpath(genpath(pwd));

%% ********************* generate the data **********************

nsample = 500;  % Number of samples,  n = 200 and 500 

p = 2000;        % Number of variables

nn = 2;         % The number of negative scale

sparsity = 10-nn;
    
v1 = zeros(p, 1);
v2 = zeros(p, 1);
v1(1:10-nn) = 1/sqrt(10)*ones(10-nn,1);
v1(11-nn:10) = -1/sqrt(10)*ones(nn,1);
v2(11:20-nn) = 1/sqrt(10)*ones(10-nn,1);
v2(21-nn:20) = -1/sqrt(10)*ones(nn,1);

randn('state',double(nsample));

P = randn(p,p);
P(:,1) = v1;
P(:,2) = v2;

P = Gram_Schmidt(P);
d = ones(p,1);
d(1:2) = [sqrt(20), sqrt(10)];
D = diag(d); 

ex_num = 200;      % Number of experiments

%% ***************** Initialization part ************************

PGD_opt1 = zeros(p,ex_num); 
PGD_opt2 = zeros(p,ex_num); 
Spnn_opt1 = zeros(p,ex_num); 
Spnn_opt2 = zeros(p,ex_num); 
PGD_time = zeros(ex_num,1); 
Spnn_time = zeros(ex_num,1); 

PGD_count1 = 0;    % record the number of correct identification of v1
PGD_count2 = 0;    % record the number of correct identification of v2
Spnn_count1 = 0;
Spnn_count2 = 0;
    
%% ************************* PGD ********************************

OPTIONS_PGD.tol = 1.0e-6;

OPTIONS_PGD.printyes = 0;

OPTIONS_PGD.maxiter = 3000;


for k = 1:ex_num
    
    randstate =(k-1)+ nsample;
    
    randn('state',double(randstate));
    
    X = randn(nsample,p);
    
    Y = X*D*P';
    
    A = Y'*Y;    % sample covariance matrix
    
   %% ********************* NPGD ***************************
    
    tstart = clock;
    
    options.tol = 1e-8;
    
    options.issym = 1;
    
    options.disp = 0;
    
    options.v0 = ones(p,1);
    
    [xint,variance] = eigs(@(y)(A*y),p,1,'LM',options);
            
    OPTIONS_PGD.Lipconst = variance;
    
    [xopt1,iter1] = PGD_PosL0_sphere(abs(xint),-A,OPTIONS_PGD,1e-5,sparsity);
    
     PGD_opt1(:,k) = xopt1;
    
    if sum(xopt1(1:10-nn)>0)==10-nn
        
        PGD_count1 = PGD_count1+1 ;   % Count the number of correct PCs
    
    end
        
    Axopt1 = A*xopt1;
    
    temp_const = xopt1'*Axopt1; 
    
    temp_A = xopt1*Axopt1';
    
    A1 = A - (temp_A +temp_A') + (temp_const*xopt1)*xopt1';
    
    options.tol = 1e-8;
    
    options.issym = 1;
    
    options.disp = 0;
    
    options.v0 = ones(p,1);
    
    [xint,variance] = eigs(@(y)(A1*y),p,1,'LM',options);
            
    OPTIONS_PGD.Lipconst = variance;
    
    [xopt2,iter2] = PGD_PosL0_sphere(abs(xint),-A1,OPTIONS_PGD,1e-5,sparsity);
    
    PGD_time(k)= etime(clock,tstart);

    if sum(xopt2(11:20-nn)>0)==10-nn
     
        PGD_count2 = PGD_count2+1 ;  %Count the number of correct PCs
    end
    
    PGD_opt2(:,k) = xopt2;
     
   %% ********************** Spnn ******************************
    
    tstart = clock;
    
    sopt1 = Spnn_solver(Y, sparsity);
    
    Spnn_opt1(:,k) = sopt1;
     
    if sum(sopt1 (1:10-nn)>0)==10-nn
       
        Spnn_count1 = Spnn_count1+1 ;   % Count the number of correct PCs
    
    end
    
    Y = Y - (Y*sopt1)*sopt1';
           
    sopt2 = Spnn_solver(Y, sparsity);
    
    Spnn_time(k) = etime(clock,tstart);
    
    Spnn_opt2(:,k) = sopt2;
    
    if sum(sopt2(11:20-nn)>0)==10-nn
        
        Spnn_count2 = Spnn_count2+1 ;  %Count the number of correct PCs
    end
       
end

%% ************ Correct rate ****************

PGD_corr1 = PGD_count1/ex_num;

Spnn_corr1 = Spnn_count1/ex_num;

PGD_corr2 = PGD_count2/ex_num;

Spnn_corr2 = Spnn_count2/ex_num;

Ptime = sum(PGD_time)/ex_num;

Stime = sum(Spnn_time)/ex_num;

%% ************************* Angle *****************************

PGD_angle1 = mean(acos(PGD_opt1'*v1))/pi*180;

Spnn_angle1 = mean(acos(Spnn_opt1'*v1))/pi*180;

PGD_angle2 = mean(acos(PGD_opt2'*v2))/pi*180;

Spnn_angle2 = mean(acos(Spnn_opt2'*v2))/pi*180;



