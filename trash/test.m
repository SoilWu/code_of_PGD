addpath('C:\Users\mayn\Documents\MATLAB\NL0_sphere\solver')
addpath('C:\Users\mayn\Documents\MATLAB\NL0_sphere\data')
addpath('C:\Users\mayn\Documents\MATLAB\lowrankquadrmax-master\matlab\bin')

%% *********** read the file *****************************

file_path = 'C:\Users\mayn\Documents\MATLAB\NL0_sphere\data\face\';

all_file_data = zeros(2429, 361);
pgm_path_list = dir(strcat(file_path, '*.pgm')); % read all file in file_path 
img_num = length(pgm_path_list); % calculate the number of file

for i = 1:img_num
    file_name = strcat(file_path, pgm_path_list(i).name); % generate the name of each file
    file_data = reshape(imread(file_name), 1, []); % read the file and reshape it
    all_file_data(i,:) = file_data; 
end

samplemean = mean(all_file_data, 1);
data = all_file_data - repmat(samplemean, 2429, 1);
A = data' * data /2428;

%% ******************** plot *********************************

I = eye(size(A,1));

figure(1)

[xopt, new_x, iter1] = PGM_solver(A, 40);
subplot(1,6,1);
imshow(mat2gray(new_x))

A = A - (xopt'*A *xopt) * (xopt * xopt');
[xopt, new_x, iter2] = PGM_solver(A, 40);
subplot(1,6,2);
imshow(mat2gray(new_x))

A = A - (xopt'*A *xopt) * (xopt * xopt');
[xopt, new_x, iter3] = PGM_solver(A, 40);
subplot(1,6,3);
imshow(mat2gray(new_x))

A = A - (xopt'*A *xopt) * (xopt * xopt');
[xopt, new_x, iter4] = PGM_solver(A, 40);
subplot(1,6,4);
imshow(mat2gray(new_x))

A = A - (xopt'*A *xopt) * (xopt * xopt');
[xopt, new_x, iter5] = PGM_solver(A, 40);
subplot(1,6,5);
imshow(mat2gray(new_x))

A = A - (xopt'*A *xopt) * (xopt * xopt');
[xopt, new_x, iter6] = PGM_solver(A, 40);
subplot(1,6,6);
imshow(mat2gray(new_x))

% figure(2)
% 
% X = all_file_data;
% 
% [xopt, new_x] = spann_solver(X, 40);
% subplot(1,6,1);
% imshow(mat2gray(new_x))


