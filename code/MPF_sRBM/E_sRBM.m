% Energy function for semi-restricted Boltzmann Machine. 
% Required for Gibbs sampling. 

% Author: Urs Koster, Jascha Sohl-Dickstein (2012)
% Web: http://redwood.berkeley.edu/wiki/Jascha_Sohl-Dickstein
% This software is made available under the Creative Commons
% Attribution-Noncommercial License.
% (http://creativecommons.org/licenses/by-nc/3.0/)

function E = E_sRBM( JW, X )
    
    % unpack weight matrices
    [nvis,T] = size(X);
    nhid = (length(JW)-nvis^2) / (nvis+1);
    J = reshape(JW(1:nvis^2), nvis, nvis);
    W = reshape(JW(nvis^2+1:end), nhid, nvis+1);
    
    % precompute common terms
    X1 = [X; ones( 1, T )]; % visible unit for hidden bias
    J = (J + J')/2; % make sure it's symmetric
    y=J*X; %ising
    ez=exp(-W*X1);  % for RBM energy term
     
    % compute energy terms, ising captures visible bias, rbm hidden bias
    E_isi = sum(X.*y); % ising model and visible bias
%%    E_rbm = sum( log(1+ez) ); % sum over visibles and samples
    E_rbm = -sum( log(1+ez) ); % sum over visibles and samples
    E = (E_isi + E_rbm) ; % energy for each data point