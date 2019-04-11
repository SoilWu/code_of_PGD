%% *********** the projection operator onto L0sphere ************
%  filename: L0sphere  
%
%% ***************************************************************
%  This code aims to seek a solution to the nonconvex nonsmooth problem  
%  
%   min 0.5||x - g||^2  s.t. ||x||_0<=kappa, ||x||=1
%    
%% ****************************************************************
  
%% ***************************************************************
%% Copyright by Wuyu Qia and Shaohua Pan, 2018/12/23
%% ***************************************************************

function xp = L0sphere(g,k)

n = size(g,1);

xsol = zeros(n,1);

xp = zeros(n,1);

gabs = abs(g);

[sg,idx] = sort(gabs,'descend'); 

sign_g = sign(g);

xsol(1:k) = sg(1:k);
    
xp(idx) = xsol/norm(xsol);

xp = xp.*sign_g;
  




