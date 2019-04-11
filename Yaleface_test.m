addpath(genpath(pwd))

%% ******************* read the file *****************************

load('Yale_face_data.mat', 'man1');   % Sample matrix

[m,n,p] = size(man1);

X = [];

for i = 1:p
    
    x = man1(:,:,i);
      
    X(:,i) = x(:);
end

X_mean = mean(X, 2);

tempX = X - repmat(X_mean,1,165);

A = tempX*tempX'/164;    % Sample covariance matrix

[nr,nc] = size(A);

%% ************ Calculate the first five PCs by PGD and SPnn **************

sparsity_list = [1000   2000  3000];

figure()

for nr = 1:3
    
    sparsity = sparsity_list(nr);
    
    options.tol = 1e-8;
    
    options.issym = 1;
    
    options.disp = 0;
    
    options.v0 = ones(nc,1);
    
    OPTIONS_PGD.tol = 1.0e-6;
    
    OPTIONS_PGD.printyes = 0;
    
    OPTIONS_PGD.maxiter = 3000;
    
    tempX = X';
    
    tempA = A;
    
    for i = 1:5
        
        [xint,variance] = eigs(@(y)(tempA*y),nc,1,'LM',options);
        
        OPTIONS_PGD.Lipconst = variance;
        
        xopt1 = PGD_PosL0_sphere(abs(xint),-tempA,OPTIONS_PGD,1e-5,sparsity);
        
        new_x1 = reshape(xopt1,m,n);
        
        %*************************** the ith image ****************************
        
        subplot(6,5,(nr-1)*10+i);
        
        imshow(mat2gray(new_x1));
        
        hold on;
        
        Axopt1 = tempA*xopt1;
        
        temp_const1 = xopt1'*Axopt1;
        
        temp_Axopt1 = xopt1*Axopt1';
        
        tempA = tempA - (temp_Axopt1 + temp_Axopt1')+ (temp_const1*xopt1)*xopt1';
        
    end
        
    % *****************calculate the first five PCs by Spnn ********************
       
    for i=1:5
        
        xopt1 = Spnn_solver(tempX, sparsity);
        
        new_x1 = reshape(xopt1, m, n);
        
        %*************************** the (i+5)th image ****************************
        
        subplot(6,5,(nr-1)*10+(i+5));
        
        imshow(mat2gray(new_x1))
        
        hold on;
        
        tempX = tempX- (tempX*xopt1)* xopt1';
    end
    
end

