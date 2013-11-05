
% try to make a general function for all models:

function high_prob = estim_highprob(p,pat_down,dat_down);

    % experiment: non-data pattern with the highest probability. 
    pn=min(5000, length(p)); % number of patterns to use
    [foo, bar] = sort(p); % value, index
    val=foo(end-pn+1:end); ind=bar(end-pn+1:end);

    % and throw out patterns that exist IRL: Create all non-data patterns:
    existo=zeros(1,pn);
    for i=1:pn % only check the biggest patterns
        if ~mod(i,pn), fprintf('.'); end
        existo(i) = sum(pat_down(ind(i))==dat_down); % does it exist?
    end
    high_prob=val(find(existo==0)); % these DONT exist in the data

return