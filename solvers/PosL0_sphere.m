%% *********** the projection operator onto PosL0+sphere ************
%  filename: PosL0_sphere  
%
%% ***************************************************************
%  This code aims to seek a solution to the nonconvex nonsmooth problem  
%  
%   min 0.5||x - g||^2  s.t. x>=0, ||x||_0<=kappa, ||x||=1
%    
%% ****************************************************************
  
%% ***************************************************************
%% Copyright by Yuqia Wu and Shaohua Pan, 2018/12/23
%% ***************************************************************

function xp = PosL0_sphere(g,k)

n = size(g,1);

xsol = zeros(n,1);

[sg,idx] = sort(g,'descend'); 

xp = [1;zeros(n-1,1)];        %%  the case sg(1)<0 is considered

if sg(k)>= 0
    
    xsol(1:k) = sg(1:k);
    
    xp(idx) = xsol /norm(xsol);
    
elseif sg(1)== 0

    abs_sg = abs(sg);
    
    u(abs_sg>0)=0;

    kk = length(find(abs_sg==0));
    
    u(abs_sg==0)=[1;zeros(kk-1,1)];
       
    xp(idx) = u;

elseif (sg(1)>0) && (sg(k)<0)
    
    indx = find(sg>=0);
    
    j = indx(end);  
    
    sgj = sg(1:j);
    
    temp_sg = [sgj; zeros(n-j,1)];
     
    xp(idx) = temp_sg/norm(sgj);     
    
end


