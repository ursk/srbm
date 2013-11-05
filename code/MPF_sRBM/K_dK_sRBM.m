% MPF objective function and gradient for semi-restricted Boltzmann Machine. , 
% Transitions allowed between all states differing by a single bit flip. 

% Author: Urs Koster, Jascha Sohl-Dickstein (2012)
% Web: http://redwood.berkeley.edu/wiki/Jascha_Sohl-Dickstein
% This software is made available under the Creative Commons
% Attribution-Noncommercial License.
% (http://creativecommons.org/licenses/by-nc/3.0/)


function [K, dK] = K_dK_sRBM( JW, X, lambda )
    
    % unpack weight matrices
    [nvis,T] = size(X);
    nhid = (length(JW)-nvis^2) / (nvis+1);
    J = reshape(JW(1:nvis^2), nvis, nvis);
    W = reshape(JW(nvis^2+1:end), nhid, nvis+1);
    
    
    % precompute common terms
    X1 = ones(nvis+1,T); X1(1:end-1,:) = X; % add always on visible unit for hidden bias
    J = (J + J')/2; % Ensure symmetric visible connections
    y=J*X;      % Ising feedforward 
    z=W*X1;     % RBM feedforward
    b=2*X-1;    % bit flips
    b1=2*X1-1;  % bit flips
    ez=exp(z);  % for RBM energy term
    
    
    % compute ising contribution to energy differences for K
    dE_isi  = 2*(b.*y  -  1/2*diag(J)*ones(1,T));
    dE_rbm = zeros( nvis, T );
    
    for n=1:nvis %
       %eW=exp(-W(:,n));
       Ai = log(1+ez);
       Bi = log( ez + exp( W(:,n)*b(n,:) ) ); % urs new derivation
       
       %%dE_rbm_n = sum(1*Ai-1*Bi)  ; % sum over hidden
       dE_rbm_n = -sum(1*Ai-1*Bi)  ; % sum over hidden
       dE_rbm(n,:) = dE_rbm_n; % rbm term for one n, all t        
    end
    
    K_nj = exp(1/2 * (dE_isi+dE_rbm) ); % terms of the MFP objective function
    K = sum(sum(K_nj)) / T;
    
    
    % compute gradient updates:
    % ising model part
    dJ =   (K_nj.*b) * X' - (1/2)*diag( sum(K_nj, 2) ) ;  % taken from previous
    dJ = (dJ + dJ')/2;
    dKdJ =  dJ / T; % update
    
    % RBM part
    dKdW = zeros( nhid, nvis+1 ); 
    dWB2n = zeros( nhid, nvis+1 );
    oh=ones(nhid,1);
    for n=1:nvis % not including bias
       %eW=exp(-W(:,n));
       dWAn = -(oh*K_nj(n,:) .* ez./(1+ez)) * X1' ; % sum over j samples happens here 
       dWB1n = ez./(ez+exp(W(:,n)*b1(n,:))) * (K_nj(n,:)'*ones(1,nvis+1) .* X1')  ; % 
       ddd =  1./(1 + exp(z - W(:,n)*b(n,:) )); % 2x100, hid*t
       dWB2n(:,n) = ddd * (K_nj(n,:).*b(n,:))'; % individual update - no bias
       dKdW = dKdW + (1*dWAn + 1*dWB1n);            % common update - bias
    end
    dKdW = dKdW + 1*dWB2n; % for some reason this term does not add up 
%%    dKdW = -1/2 * dKdW / T; % 1/2 is pulled down from K=exp(1/2 ediff)
    dKdW = 1/2 * dKdW / T; % 1/2 is pulled down from K=exp(1/2 ediff)
    
    %concatenate and vectorize all gradients
    dK = [dKdJ(:); dKdW(:)]; % 
    
    K = K + lambda * sum(abs(JW(:))); % L1 regularizer
    dK = dK + lambda * sign( JW(:) );
    %K = K + lambda * sum(JW(:).^2); % L2 regularizer
    %dK = dK + 2*lambda * ( JW(:) );
return
