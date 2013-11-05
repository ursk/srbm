function [Wrbm, Z_rbm] = estim_rbm(X, N, d_vis, d_hid, minf_options, lambda, PMPF, Winit)
    
    % parameter estimation
    
    fprintf('RBM model\n')
    %if nargin==7
        Winit = 0.1 * randn( d_hid+1, d_vis+1 ) / sqrt(d_vis+1);
        Winit(end,1:end-1)=2;
    %end
    t_min=tic();
    if PMPF
        minf_options.LINESEARCH_NOTIFY = 1;           
        sample_iter = 10; % Gibbs steps -- try less steps to speed up. 
        sample_n_linesearches = 1; % no subsampling but full gibbs
        include_standard_MPF = 0; % actually much worse with this, and it's only.   
        Wmpfn = minFunc_fs( @K_dK_RBM_PMPF_subsample, Winit(:), minf_options, X, lambda, sample_iter, sample_n_linesearches, include_standard_MPF ); %% estimation (expensive)
    else
        Wmpfn = minFunc( @K_dK_RBM, Winit(:), minf_options, X, lambda ); %% estimation (expensive)
    end
    t_min = toc(t_min);
    fprintf( 'RMB parameter estimation in %f seconds \n', t_min );
    Wrbm = reshape( Wmpfn, size( Winit ) );
    
    % tail matching
    [best_state_vectors, best_state_probs, rndpr] = get_best_states(X, N);
    T=size(X,2);
    for i=1:N
        %best_state_probs(i) = sum(downstate == rndpr*best_state_vectors(:,i)) / T;
        xx = best_state_vectors(:,i);
        %p0true = sum(rndpr* X==rndpr*xx) / T; % use ground truth to calibrate
        p0raw_rbm = exp(-E_RBM(Wrbm, xx ));
        Z_rbm(i) = p0raw_rbm / best_state_probs(i);
    end
    Z_rbm  = sum(Z_rbm .* best_state_probs) / sum(best_state_probs);
    
return
