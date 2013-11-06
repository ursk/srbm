function [x_out, logZ_e, Z, logaisS, logZ_all] = ais_converge_wrapper(rate, J, W, JW, p, desc, plot_ind)

logZ_all = [];

logaisS = 0:ceil(log2(p.aisS));

if nargin == 6  
    plot_ind = ceil( rand()*100000);
end

%rate(:) = 0.5; % start at a uniform distribution rather than independent model  

for ind = 1:length( logaisS )
    p.aisS = ceil(2^( logaisS(ind) ));

    [x_out, logZ_e, Z] = ais(rate, J, W, JW, p); % Annealed Importance Sampling
    logZ_all = [logZ_all, logZ_e(:)];
end    
% moved outside loop because it leaks files and causes too many files open    
    hand = figure(plot_ind); clf;
    semilogx( 2.^(logaisS(1:ind))', log2(exp(logZ_all)), '.' ); hold on;
    semilogx( 2.^(logaisS(1:ind)), log2(mean(exp(logZ_all))) ); % should use log2 here for bits!!
    xlabel( 'Number of intermediate distributions' );
    ylabel( 'log2 Z estimate' );
    try
        ylim([log2(mean(exp(logZ_e(:))))-1/50   log2(mean(exp(logZ_e(:))))+1/50]), %grid on %% 1/50 is an error of 1 bit
    catch me
        me
    end
    title( sprintf( 'Convergence of AIS for %s', desc) );
    drawnow;
    try
        saveas( hand, sprintf( 'models/ais/%s.pdf', desc) );
        % saveas( hand, sprintf( 'models/ais/%s.fig', desc) );
        % can't save figures so save as mat ok pls
        % save(sprintf( 'models/ais/%s.mat', desc))
    catch me
        me
    end
    




