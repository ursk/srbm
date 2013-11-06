% Demo script to estimate all three model types (Ising, RBM, sRBM) using MFP
% and evaluate model probabilities using AIS

function bits = cell_estim()
    ticall=tic();
    
    addpath('MPF_ising')
    addpath('MPF_RBM')
    addpath('MPF_sRBM')
    addpath(genpath('3rd_party_code/minFunc/minFunc'))

    % --------------------------
    % set model parameters
    % --------------------------
    
    modelname='ising_rbm_srbm_demo';
    lambda = 0.001;
    seed = 42;
    timesteps=1;
        
    SPEEDUP=5; % if >1, train a smaller model (faster, less accurate) for debugging
    dataseed = seed; % randomize everything again. 
    PMPF_for_RBM = 0 % new implementation Jascha made -- we don't really want this anymore. 
    
    p.t      = 1;    % number temporal units
    p.bin    = 3;    % binning of consecutive 6.6ms frames
    p.d_hid  = 25;   % number hidden units 
    p.d_s    = 25;   % number hidden units for sRBM
    p.T      = 50e3/SPEEDUP; % size of training data
    p.lambda = lambda;      % regularization parameter for L1
    plot_ind = 100 * ceil( rand()*1000); % figure number for AIS convergence plots.
    p.burnin = 10; % sampling steps for gibbs. 
    
    p.s      = 21; % number of spatial units   
    p.d      = p.s*p.t;  % total model size
    p.aisT   = 500;   % number of AIS samples
    p.aisS   = 1e5/SPEEDUP; % number of AIS annealing steps 
    
    minf_options = [];
    minf_options.maxFunEvals = 5000; % use a lot, just cause. Still 500 Iterations! 
    minf_options.maxIter = 1000/SPEEDUP; % [MODIFIED] should be 200 
    minf_options.progTol = 1e-12; % hit defaule 1e-9 too easily...
    minf_options.optTol = 1e-12; % hit defaule 1e-9 too easily...
    minf_options.TolX = 1e-12; % 1e-9 maybe important for sRBM
    minf_options.Display = 'iter'; % default is 'iter' or 'final'
    
    minf_options2=minf_options;
    minf_options2.maxIter = minf_options.maxIter/2;
    
    [foo bar]=system('hostname');
    exactZ    = (p.d<=23) & strcmp(strtrim(bar),'htpc2') % we can compute the exact Z for up to 23 units
    
    % --------------------------
    % create a dummy data set
    % --------------------------
    % for demonstration, create a random data set with higher order dependencies.
    rand('seed',dataseed), randn('seed',dataseed) % one fixed dataset, randomize only parameters and AIS
    load('demodata.mat', 'XY')
    rand('seed',seed), randn('seed',seed) % different seed for AIS and parametres
    X = XY(:,1:p.T); % training
    Y = XY(:,p.T+1:2*p.T); % validation
    
  
    rndprj=rand(1,p.d);
    dat_down=rndprj*Y; % used for p_emp and scatter plot
    
    % estimate the number of states for Z-matching
    [dummy, order] = sort(rndprj*X); % do projection on X, not Y, for this! 
    transitions    = find(diff(dummy));
    p.N            = 10; % size(transitions,2) - 20; % tail matching size, can be more than we have unique states. some get lost?
    fprintf('Tail matching with N=%d-20 terms in the expansion\n', p.N+20)
    
    if exactZ
        fprintf('Enumerating all %d patterns for exact Z ...', 2^p.d)
        tic
        patterns = all_states(p.d);
        pat_down=rndprj*patterns';
        fprintf('took %d seconds', toc)
    end
    
    
    % --------------------------
    % estimate models 
    % --------------------------    
    
    p
    fprintf(['\n\n  ESTIMATING MODEL ', modelname, ' \n\n'])
    
    % empirical
    fprintf('Counting empirical probabilities on %d samples...\n', length(X))
    pemp = estim_empirical(dat_down); % empirical distribution
    
    
    % independent model
    [rate, Z_ind_TM7] = estim_ind(X, p.N);
    rate = rate + 1/p.T; % [HACK] make sure no rates are zero so we don't get -inf log prob
    pind = prod( bsxfun(@times, rate, Y) + bsxfun(@times,(1-rate),(1-Y))  ); % product of normalized probs
    
    
    
    % 
    for iii = [1 2 3]  % 2 3   % cell_estim('aisplot', 's', 0.001, 0)     !head -n 200 cell_estim.m
    
    if iii==1;   
        %ISING
        fprintf('Estimating ISING MODEL')
        [Jisi, Z_tm_isi, Z_i, p_i] = estim_isi(X, p.N, minf_options, p.lambda*10/10, 'allflip'); %style ='allflip'; % Needs stronger regularization
        [x_out, logZ_e_isi, Z_sam_isi, logaisS, logZ_all_isi] = ais_converge_wrapper(rate, Jisi, [], [], p, [modelname,'_isi'],plot_ind+1); % Annealed Importance Sampling
        %x_out=0, logZ_e_isi=0, Z_sam_isi=0
        pisi = 1/Z_sam_isi * exp(-sum( (Jisi*Y).*Y, 1 ));
        mi_isi = mean(log2(pisi) - log2(pind));
        
        fprintf('----------------------\n')
        fprintf('All-flip Ising has %2.2f bits/s\n', 150/p.bin*mi_isi); 
        bits(iii)=150/p.bin*mi_isi; 
        fprintf('----------------------\n')
        if exactZ
            p_mo = exp(-sum( (Jisi*patterns').*patterns', 1 ))/ 1; % model probabilities. Was /Z_isi, but want raw Z.
            Z_emp_isi = sum (p_mo);
            pisi = 1/Z_emp_isi * exp(-sum( (Jisi*Y).*Y, 1 ));
            mi_isi_true = mean(log2(pisi) - log2(pind));
            high_prob_isi = estim_highprob(p_mo,pat_down,dat_down);
            fprintf('\nCOMPARISON OF PARTITION FUNCTION:\n')
            fprintf('Z_exact = %2.3f\n', Z_emp_isi)
            fprintf('Z_sampling = %2.3f\n', Z_sam_isi)
            fprintf('Z_tm = %2.3f\n', Z_tm_isi)
            fprintf('--------------------------------\n')
        else
        end


    end  
    
    if iii==2;
        % RBM
        %convolution = 1; 
        fprintf('Estimating RBM MODEL') 
        %minf_options.maxIter = minf_options.maxIter/2;
        
        % -- first half of RBM -
        Winit = 0.001 * randn( p.d_hid+1, p.d+1 ) / sqrt(p.d+1);
        [Wfirst, Z_tm_rbm] = estim_rbm(X,  p.N, p.d, p.d_hid, minf_options2, p.lambda, PMPF_for_RBM, Winit); % no convolution
        % --
       
        [Wrbm, Z_tm_rbm] = estim_rbm(X,  p.N, p.d, p.d_hid, minf_options2, p.lambda, PMPF_for_RBM, Wfirst); % two step estimation!
        [x_out, logZ_e_rbm, Z_sam_rbm, logaisS, logZ_all_rbm] = ais_converge_wrapper(rate, [], Wrbm, [], p, [modelname,'_rbm'],plot_ind+2); % Annealed Importance Sampling
        prbm = 1/Z_sam_rbm * exp(-E_RBM(Wrbm, Y ));
        mi_rbm =  mean(log2(prbm) -log2(pind)) ;
        fprintf('----------------------\n')
        fprintf('RBM has %2.2f bits/s\n', 150/p.bin*mi_rbm); 
        bits(iii)=150/p.bin*mi_rbm; 
        fprintf('----------------------\n')
        if exactZ
            p_mo = exp(-E_RBM(Wrbm, patterns' )) ;
            Z_emp_rbm = sum(p_mo) ;
            prbm = 1/Z_emp_rbm * exp(-E_RBM(Wrbm, Y ));
            mi_rbm_true =  mean(log2(prbm) -log2(pind)) ;
            high_prob_rbm = estim_highprob(p_mo,pat_down,dat_down);
            fprintf('\nCOMPARISON OF PARTITION FUNCTION:\n')
            fprintf('Z_exact = %2.2f\n', Z_emp_rbm)
            fprintf('Z_sampling = %2.2f\n', Z_sam_rbm)
            fprintf('Z_tm = %2.2f\n', Z_tm_rbm)
            fprintf('--------------------------------\n')
        else
        end
    end
    
    if iii==3
        % SRBM
        fprintf('Estimating sRBM MODEL using only RBM initialization...') 
        bias_init =  diag(Wfirst(p.d_s+1,1:p.d));
        W_init    = Wfirst(1:p.d_s,:);
        [JWout, Z_tm_srbm] = estim_srbm(X, p.N, p.d, p.d_s, minf_options2, p.lambda, bias_init, W_init); % strip off vis bias -- 
        [x_out, logZ_e_srbm, Z_sam_srbm, logaisS, logZ_all_srbm] = ais_converge_wrapper(rate, [], [], JWout, p, [modelname,'_srbm'],plot_ind+3);
        psrbm = 1/Z_sam_srbm * exp(-E_sRBM(JWout, Y ));
        mi_srbm = mean(log2(psrbm) - log2(pind)) ; 
        fprintf('----------------------\n')
        fprintf('sRBM has %2.2f bits/s\n', 150/p.bin*mi_srbm);
        bits(iii)=150/p.bin*mi_srbm; 
        fprintf('----------------------\n')
        if exactZ
            p_mo = exp(-E_sRBM(JWout, patterns' ));
            Z_emp_srbm = sum(p_mo);
            psrbm = 1/Z_emp_srbm * exp(-E_sRBM(JWout, Y ));
            mi_srbm_true = mean(log2(psrbm) - log2(pind)) ;
            high_prob_srbm = estim_highprob(p_mo,pat_down,dat_down);
            fprintf('\nCOMPARISON OF PARTITION FUNCTION:\n')
            fprintf('Z_exact = %2.2f\n', Z_emp_srbm)
            fprintf('Z_sampling = %2.2f\n', Z_sam_srbm)
            fprintf('Z_tm = %2.2f\n', Z_tm_srbm)
            fprintf('--------------------------------\n')
        else
        end
    end
    end % parfor
    
    
    if exactZ
        % convergence check for AIS 
        clf
        subplot(2,2,1)
        if exist('logZ_all_isi')    
            plot(log2(exp(logZ_all_isi))','.'), hold on
            plot(log2(mean(exp(logZ_all_isi)))','-r','LineWidth', 2), hold off
            line([0 20], [log2(Z_emp_isi) log2(Z_emp_isi)])
            xlim([0.5 18.5]), ylim([log2(Z_emp_isi)-.1 log2(Z_emp_isi)+.1])
        end
        subplot(2,2,2)
        if exist('logZ_all_rbm')   
            plot(log2(exp(logZ_all_rbm))','.'), hold on
            plot(log2(mean(exp(logZ_all_rbm)))','-r','LineWidth', 2), hold off
            line([0 20], [log2(Z_emp_rbm) log2(Z_emp_rbm)])
            xlim([0.5 18.5]), ylim([log2(Z_emp_rbm)-.1 log2(Z_emp_rbm)+.1])
        end
        subplot(2,2,3)
        if exist('logZ_all_srbm')
            plot(log2(exp(logZ_all_srbm))','.'), hold on
            plot(log2(mean(exp(logZ_all_srbm)))','-r','LineWidth', 2), hold off
            line([0 20], [log2(Z_emp_srbm) log2(Z_emp_srbm)])
            xlim([0.5 18.5]), ylim([log2(Z_emp_srbm)-.1 log2(Z_emp_srbm)+.1])
        end
        subplot(2,2,4)
        if exist('logZ_all_gibbs')
            plot(log2(exp(logZ_all_gibbs))','.'), hold on 
            plot(log2(mean(exp(logZ_all_gibbs)))','-r','LineWidth', 2), hold off
            line([0 20], [log2(Z_emp_gibbs) log2(Z_emp_gibbs)])    
            xlim([0.5 18.5]), ylim([log2(Z_emp_gibbs)-.1 log2(Z_emp_gibbs)+.1])
        end
        saveas(gcf, ['./exactZ/exactZcomparison-', modelname], 'pdf')   
    else
    end
    
    
    
    % save results 
    clear p_mo pat_down X XY Xall Xb Xbin Xsub Y patterns % delete large variables
    fprintf('Writing output to %s \n', modelname);
    save(['./', modelname, '.mat']) 
    fprintf('Total runtim was %d minutes \n', round(toc(ticall)/60))
return

