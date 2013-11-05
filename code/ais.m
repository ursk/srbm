% Annealed Importance Sampling.
% [DONE] AIS for Ising model, -- Ice ice baby!
% [TODO]  generatize to RBM and sRBM 
%
    % Using pseudocode from Jascha
    % FOR n_samples
    %   x = draw new sample from independent model
    %   logZ_es = log Z from independent model
    %   FOR beta annealing steps    
    %       logZes += [mixed energy with current beta at current x] - [mixed energy with previous beta at current x]
    %       SHOULD BE: logZes += old - new
    %       sample new x using [mixed energy with current beta]
    %   end
    %   store logZ_es into logZ_e
    % end
    % Z = mean exp (logZ_e)
    
    
function [x_out, logZ_e, Z] = ais(rate, J, W, JW, p)
    tic
    %Indep=inline('-sum(log(  repmat(rate,1,size(x,2)).*x +  repmat(1-rate, 1,size(x,2)).*(1-x)  ))'); % log-probablity (normalized)
    %Ising =inline('sum(X .* (J*X))'); % energy
    %sigmoid=inline('1./(1+exp(-x))');
    
    T = p.aisT;    % samples to draw
    s = p.aisS;    % steps for each sample
    d = length(rate); % data dimensions
    %uni_rnd=ceil(d*rand(T*s,1)); % precompute random -- THIS BECOMES HUUUUUUUUUGE! (and it's never used)
    xout=zeros(d,T);
    nsamples=1; % for SW
    SW=0; % don't actually use SW, get it to work with Gibbs first
    
    x = double(rand(d,T) < repmat(rate,1,T)); % T samples in parallel
    logZ_es = log(1); % the indep model comes normalized, this is 0
    
    % --------------------
    %%% Big, slow loop %%%
    % --------------------
    for b = [1:s]/s % annealing steps, one bit flip per step
        if ~mod(b*s,10000); fprintf('(%dk)', s*b/1000);end
        % AIS step
        if length(J) % ising
            E_m = Ising(J, x); %[!!!] probably not really redundant with bit flip since all 100 x are computed at once. 
        elseif length(W) % RBM
            theta = W(:); 
            E_m = E_RBM(theta, x); % 
        else % sRBM
            E_m = E_sRBM(JW, x);
        end
        E_i = Indep(rate, x); %[!!!] would need clever indexing, not worth it for a 33% speedup. 
        E_old = (b-1/s)*E_m + (1-b+1/s)*E_i;  % at old b -- CURRENT X.
        E_new = b*E_m + (1-b)*E_i; % mixed energy at current b
        
        % the imporatance weights w=f1/f2... update step, in log world
        logZ_es = logZ_es + (-E_new + E_old); % 5000 energies -> 5000 logZ's
        
        % GIBBS: sample new x (at new mixed E)
        n=ceil(d*rand(T,1)); % randomly choose dimension to update
        iae= (d*[1:T]-d)' + n; % index active elements in all 100 vectors
        x_test0=x;      x_test1=x; 
        x_test0(iae)=0; x_test1(iae)=1;                
        
        if length(J)
            E_m0 = Ising(J, x_test0); % [!!!]
            E_m1 = Ising(J, x_test1); % [!!!]
        elseif length(W)
            theta = W(:); 
            E_m0 = E_RBM(theta, x_test0);
            E_m1 = E_RBM(theta, x_test1);
        else % sRBM
            E_m0 = E_sRBM(JW, x_test0);
            E_m1 = E_sRBM(JW, x_test1);
        end
        E0 = b*E_m0+(1-b)*Indep(rate, x_test0); % [!!!]
        E1 = b*E_m1+(1-b)*Indep(rate, x_test1); % [!!!]
        
        %xvec = rand(1,T) < sigmoid(-E1 + E0); % samples under transition kernel
        xvec = rand(1,T) < 1./(1+exp(-(-E1 + E0)));
        newstate=zeros(1,T); % default:set to zero
        newstate(find(xvec)) = 1; % elements where we exceed threshold: set to one
        x(iae) = newstate; % set the correct elements to the new state
        
    end
    logZ_e = logZ_es; % these are all in one now, no more k 
    x_out =x;
    Z = mean( exp(logZ_e) );
    fprintf('\nSampling %d times %d steps took %2.2fs\n', T, s, toc)

return

% dont use inline() since it seems really slow. 

% Indep is really expensive!  
function out = Indep(rate, x) % x is 25x100, i.e. cells x samples 
    a = bsxfun(@times,rate,x) +  bsxfun(@times,(1-rate),(1-x)); % this is expensive, 100k calls, WHY???
    out = -sum(log(  a  )); % log-probablity (normalized)

function out = Ising(J,X) 
    out = sum(X .* (J*X)); % energy, this is now also becoming expensive. 

%function out = sigmoid(x)
%    out =  1./(1+exp(-x));









    
    % loop independent samples 
    % parfor k=1:T
    %     fprintf('.')
    %     x = double(rand(d,1) < rate); % sample independent model
    %     logZ_es = log(1); % the indep model comes normalized, this is 0
    % 
    %     for b = [1:s]/s % annealing steps, one bit flip per step
    %         
    %         % AIS step
    %         if length(J) % ising
    %             E_m = Ising(J, x);
    %         elseif length(W) % RBM
    %             theta = W(:); 
    %             E_m = E_RBM(theta, x); % ILLEGAL CALL -- access W(0,11)
    %         else % sRBM
    %             E_m = E_sRBM(JW, x);
    %         end
    %         E_i = Indep(rate, x);
    %         E_old = (b-1/s)*E_m + (1-b+1/s)*E_i;  % at old b -- CURRENT X.
    %         E_new = b*E_m + (1-b)*E_i; % mixed energy at current b
    %         
    %         logZ_es = logZ_es + (-E_new + E_old); % this had a sign error???
    %         
    %         % GIBBS: sample new x
    %         n=ceil(d*rand()); % randomly choose dimension to update
    %         x_test0=x; x_test1=x; 
    %         x_test0(n)=0; x_test1(n)=1;                
    %         
    %         if length(J)
    %             E_m0 = Ising(J, x_test0);
    %             E_m1 = Ising(J, x_test1);
    %         elseif length(W)
    %             theta = W(:); 
    %             E_m0 = E_RBM(theta, x_test0);
    %             E_m1 = E_RBM(theta, x_test1);
    %         else % sRBM
    %             E_m0 = E_sRBM(JW, x_test0);
    %             E_m1 = E_sRBM(JW, x_test1);
    %         end
    %         E0 = b*E_m0+(1-b)*Indep(rate, x_test0); % energy of flipped state 
    %         E1 = b*E_m1+(1-b)*Indep(rate, x_test1); % has higher energy
    %         
    %         x(n) = 0; 
    %         if rand() < sigmoid(-E1 + E0); % samples under transition kernel
    %             x(n) = 1; 
    %         end %
    %         
    %         % SWENDSEN-WANG: sample new x and do it better
    %         % if SW
    %         %     d_vis=d; d_hid=0;
    %         %     J_joint = zeros( d_vis + d_hid );
    %         %     J_joint(1:d_vis, 1:d_vis) = J; % visible to visible
    %         %     J_joint( d_vis+1:end, 1:d_vis ) = W(:,1:d_vis)/2; % visible to hidden
    %         %     J_joint( 1:d_vis, d_vis+1:end ) = W(:,1:d_vis)'/2;% visible to hidden
    %         %     J_joint = J_joint + diag( [zeros(d_vis,1); W(:,end)] ); % hidden biases
    %         % 
    %         %     X_joint = round(rand(d_vis+d_hid,nsamples));
    %         %     X_joint = sw_allall_bm( X_joint, -(J_joint - diag(diag(J_joint))), -(diag(J_joint)), sample_iter );
    %         %     X = X_joint( 1:d_vis, : );
    %         % end
    %     end
    %     logZ_e(k) = logZ_es;
    %     x_out(:,k)=x;
    % end
    %Z = mean( exp(logZ_e) );
    %fprintf('\nSampling %d times %d steps took %2.2fs\n', T, s, toc)