function patterns = all_states(d)

    fprintf('generating patterns (Jaschas method)...')
    tic
    for ii = 1:d
        patterns(:,ii) = bitget( (0:(2^d-1))', ii );
    end
    fprintf('done in %2.0fs\n', toc)
    
return