function [xopt] = sparse_nonnegative(g, k)

n = size(g,1);

x = zeros(n,1);

[sg,idx] = sort(g,'descend'); 

if sg(k)>=0
    
    x(1:k) = sg(1:k);
    
else 
  
    for i = 1: k-1
        
        if sg(k-i)>=0
            
            x(1:k-i) = sg(1:k-i);
            
            break
        end
    end
end

xopt = zeros(n, 1);

xopt(idx) = x;

xopt = xopt/norm(xopt);