% do not call this with patterns for X, it does not make sense

function [pemp] = estim_empirical(downstate)
    
    T = length(downstate);
    
    fprintf('counting %d states...', T)
    tic
    pemp=zeros(T,1);
    for i=1:T % 1M states to visit, need to perform a serach in each. 
        if ~mod(i,1000); fprintf('.'); end
        pemp(i) = sum( downstate(i)==downstate)/T;
    end
    fprintf('done in %2.0fs\n', toc)
    pemp=pemp';
return
 