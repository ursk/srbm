

function [JWout, Z_srbm] = estim_srbm(X, N, d_vis, d_hid, minf_options, lambda, Jinit, Winit)
    
    % parameter estimation
    
    if nargin==6; 
        fprintf('\n\n\n\nInitializing sRBM with random\n\n'); 
        Jinit = 0.001 * randn(d_vis); Jinit=Jinit+Jinit'; %ising
        Winit = 0.01 * randn( d_hid+0, d_vis+1 ) / sqrt(d_vis+1); %rbm
    else; 
        fprintf('\n\n\n\nInitializing sRBM with Ising and RBM\n\n'); 
        %Jinit = 0.001 * randn(d_vis); Jinit=Jinit+Jinit'; %ising
    end
    
    JWinit = [Jinit(:); Winit(:)];
    t_min=tic();
    JWout = minFunc( @K_dK_sRBM, JWinit(:), minf_options, X, lambda ); %% estimation (expensive)
    t_min = toc(t_min);
    fprintf( 'sRBM parameter estimation in %f seconds \n', t_min );
    
    Jout = reshape(JWout(1:d_vis^2), d_vis, d_vis);
    Wout = reshape(JWout(d_vis^2+1:end), d_hid+0, d_vis+1);

    % partition function estimation
    
    T=size(X,2);
    [best_state_vectors, best_state_probs, rndpr] = get_best_states(X, N);
    for i=1:N % TM_7 estimation
        xx = best_state_vectors(:,i);
        %p0true = sum(rndpr* X==rndpr*xx) / T; % use ground truth to calibrate
        p0raw_srbm = exp(-E_sRBM(JWout,xx ));
        Z_srbm(i) = p0raw_srbm / best_state_probs(i);
    end
    Z_srbm = sum(Z_srbm .* best_state_probs/sum(best_state_probs)); 
    
    
return