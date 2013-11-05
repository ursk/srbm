function [Jisi, Z, Zi, best_state_probs] = estim_isi(X, N, minf_options, lambda, allflip)
    
    % perform parameter estimation
    
    T = size(X,2);
    d_vis =size(X,1);
    Jinit = randn(d_vis) / sqrt(d_vis) / 100; % biases on the diagonal
    Jinit = (Jinit+Jinit')/2; % 
    t_min = tic();
    if ~exist('allflip')
        fprintf( '\nRunning %d iterations minFunc for Ising model...\n', minf_options.maxIter );
        J = minFunc( @K_dK_ising, Jinit(:), minf_options, X, lambda );
    elseif strcmp(allflip,'allflip')
        fprintf( '\nRunning %d iterations minFunc for allflip-Ising model...\n', minf_options.maxIter );
        J = minFunc( @K_dK_ising_allbitflipextension, Jinit(:), minf_options, X, lambda );
    end
    Jisi = reshape(J, size(Jinit));
    t_min = toc(t_min);
    fprintf( 'Ising parameter estimation in %f seconds \n', t_min );
    
    % tail matching
    
    [best_state_vectors, best_state_probs, rndpr] = get_best_states(X, N);
    for i=1:N
      xx = best_state_vectors(:,i); % how are these best if there are zeros?
      %p0true = sum(rndpr* X==rndpr*xx) / T; % use ground truth to calibrate. This is a count, can it be zero??
      p0raw_isi = exp(-sum( (Jisi*xx).*xx, 1 )); 
      Zi(i) = p0raw_isi / best_state_probs(i); % becomes infinity if p0 is zero
    end
    %fprintf('DEBUG MESSAGE Z terms are\n')
    %Zi(1:10)
    Z  = sum(Zi .* best_state_probs) / sum(best_state_probs);
     
     
return
