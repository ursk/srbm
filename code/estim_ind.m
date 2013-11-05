function [rate, Z_ind] = estim_ind(X, N)
    
    % independent model estimation
    T = size(X,2);
    rate = mean(X,2);
    
    % partition function estimation -- not needed, this model is easy
    [best_state_vectors, best_state_probs, rndpr] = get_best_states(X, N);
    for i=1:N % TM_7 estimation
        xx = best_state_vectors(:,i);
        %p0true = sum(rndpr* X==rndpr*xx) / T; % use ground truth to calibrate. would this not be best_stat_probs???
        p0raw_ind = prod( bsxfun(@times, rate, xx) + bsxfun(@times,(1-rate),(1-xx)) ); % not same as true! 
        Z_ind(i) = p0raw_ind / best_state_probs(i);
    end
    Z_ind  = sum(Z_ind .* best_state_probs/sum(best_state_probs)); 
    
return