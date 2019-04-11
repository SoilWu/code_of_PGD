addpath('C:\Users\mayn\Documents\MATLAB\NL0_sphere\solver')
addpath('C:\Users\mayn\Documents\MATLAB\NL0_sphere\data')
addpath('C:\Users\mayn\Documents\MATLAB\lowrankquadrmax-master\matlab\bin')

% for i = 1:5
%     subplot(6,5,i)
%     imshow(mat2gray(OneA(:,100*(i-1)+1:100*i)))
%     subplot(6,5,5+i)
%     imshow(mat2gray(TwoA(:,100*(i-1)+1:100*i)))
%     subplot(6,5,10+i)
%     imshow(mat2gray(ThereeA(:,100*(i-1)+1:100*i)))
%     subplot(6,5,15+i)
%     imshow(mat2gray(OneB(:,100*(i-1)+1:100*i)))
%     subplot(6,5,20+i)
%     imshow(mat2gray(TwoB(:,100*(i-1)+1:100*i)))
%     subplot(6,5,25+i)
%     imshow(mat2gray(ThreeB(:,100*(i-1)+1:100*i)))
% end



 gap = 10;
% bottle1 = 0.05 * ones(810,500+4*gap);


% for i = 1:5
%     bottle1(1:100,100*(i-1)+1+gap*(i-1):100*i+gap*(i-1)) = standard(:,100*(i-1)+1:100*i);
% end

% for i = 1:5
%     bottle1(1:100,100*(i-1)+1+gap*(i-1):100*i+gap*(i-1)) = OneA(:,100*(i-1)+1:100*i);
% end
% 
% for i = 1:5
%     bottle1(131:230,100*(i-1)+1+gap*(i-1):100*i+gap*(i-1)) = OneB(:,100*(i-1)+1:100*i);
% end
% 
% for i = 1:5
%     bottle1(291:390,100*(i-1)+1+gap*(i-1):100*i+gap*(i-1)) = TwoA(:,100*(i-1)+1:100*i);
% end
% 
% for i = 1:5
%     bottle1(421:520,100*(i-1)+1+gap*(i-1):100*i+gap*(i-1)) = TwoB(:,100*(i-1)+1:100*i);
% end
% 
% for i = 1:5
%     bottle1(581:680,100*(i-1)+1+gap*(i-1):100*i+gap*(i-1)) = ThereeA(:,100*(i-1)+1:100*i);
% end
% 
% for i = 1:5
%     bottle1(711:810,100*(i-1)+1+gap*(i-1):100*i+gap*(i-1)) = ThreeB(:,100*(i-1)+1:100*i);
% end
% 
% figure(1)
% 
% imshow(mat2gray(bottle1))



% for i = 1:5
%     V(:,100*(i-1)+1:100*i) = man1(:,:,i);
% end
% 
% 
% 
% standard = zeros(100, 500);
%  
% for i = 1:5
%     standard(:,100*(i-1)+1:100*i) = reshape(U(:,i),100,100);
% end
% 
% bottle = 0.050 *ones(220,500+4*gap);
% 
% I = 0.0002*eye(100);
% 
% for i = 1:5
%     bottle(1:100,100*(i-1)+1+gap*(i-1):100*i+gap*(i-1)) = V(:,100*(i-1)+1:100*i)*I;
% end
% 
% for i = 1:5
%     bottle(121:220, 100*(i-1)+1+gap*(i-1):100*i+gap*(i-1)) = standard(:,100*(i-1)+1:100*i);
% end
% 
% figure(1)
% imshow(mat2gray(bottle))

% [m, n, p] = size(man1);
% X = [];
% for i = 1:p
%     x = man1(:,:,i);
%     x = reshape(x, [], 1);
%     X(:, i) = x;
% end
% 
% X_mean = mean(X, 2);
% X = X - repmat(X_mean, 1, 165);
% A = X * X'/164;
% 
% v1 = U(:,1);
% variance = v1'*A*v1;
% sparsity = linspace(500,5000,10);
% 
% y1 = [];
% y2 = [];
% 
% for i = 1:10
%     [xopt, iter1] = PGM_solver(A, sparsity(i));
%     y1(i) = xopt'*A*xopt/variance;
%     
%     [xopt] = Spnn_solver(X', sparsity(i));
%     y2(i) = xopt'*A*xopt/variance;
% end
% 
% plot(sparsity,y1,'r-o',sparsity,y2,'b-')
    

