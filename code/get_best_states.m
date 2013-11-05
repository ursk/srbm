% function is called from estim_in and estim_rbm
function [best_state_vectors, best_state_probs, rndpr] = get_best_states(X, N, plotthis)
    
    T=size(X,2); % 10k
    d=size(X,1); % 20
    rndpr=rand(1,d);
    downstate =  rndpr* X; % project data onto random vector
    [srtdata, order] = sort(downstate); %10k ordered data states
    
    Xsort = X(:, order); % by random, but group identical states
    
    places = [0, find(diff(srtdata))]+1; % 1....; 3379...;  3649...; are 
    [sizes, sizeorder] = sort(-diff(find([1,diff(srtdata)]))); % each place where a transition occurs, how long is it?
    
    fprintf('TAILMATCHING: Trying to extract %d states from %d total\n', N, length(sizeorder))    
    best_state_vectors = Xsort(:, places(sizeorder(1:N))); % this fails?
    best_state_probs = -sizes(1:N) / T;
    
    
    
    % transitions = find(diff(srtdata)); % index of transitions in sorted space
    % [sorted, indexed] = sort(diff(transitions)); % index where biggest transitions occur, biggest at the end
    % 
    % % plot power law distribution of cells
    % if exist('plotthis')
    %     figure(99)
    %     loglog(fliplr(sorted),'.-'), xlabel('state'), ylabel('count') % now we see the powah law!! 
    %     line([300 1], [1 15000]) % fit a line -- slope m=50 wow that's nasty! but maybe m>10 is a candidate for this method. 
    %     title('figure 2 - power law staes')
    % end
    % 
    % % identify the big states, then use them, weighted by probability, to compute Z.
    % Xsort = X(:, order); % by random, but group identical states
    % size(indexed)
    % best=transitions(indexed(end-N+2:end)); % index of where a large jump occurs in ordered data
    % counts = sorted(indexed(end-N+2:end)); % sawtooth pattern, wtf?
    % best_state_vectors = Xsort(:, [best+1, 1]); % last vector is zeros state. Can contain states that occor 0 times???
    % 
    % for i=1:N
    %     best_state_probs(i) = sum(downstate == rndpr*best_state_vectors(:,i)) / T;
    % end
    % 
    % keyboard
    
    
    
    
return
    