
function [pxn] = decodeBin(time_bin_start, time_bin_end, spiketimes, fr_per_bin, opts)

    % Get number of spikes
    n = cellfun(@(x) sum(x >= time_bin_start & x < time_bin_end), spiketimes);
    
    % Check number of active cells is at least 4
    if sum(n > 0) >= 4
        
        % Calculate P(n|x), log-probability of spikes assuming Poisson distribution for each binned position
        for i = 1:length(spiketimes)
            pnix(:,i) = log(pdf('Poisson', n(i), fr_per_bin(:,i) * opts.binsize));
        end
        pnix(isinf(pnix)) = -100; % Cap low probabilities
        pnx = sum(pnix,2);
        
        % Calculate P(x), log-probability of every position; assume uniform
        px = log(ones(1,size(fr_per_bin,1))/size(fr_per_bin,1));
        
        % Calculate P(n), log-probability of observed spike vector; assume Poisson firing rates
        for i = 1:length(spiketimes)
            pni(i) = log(pdf('Poisson', n(i), nanmean(fr_per_bin(:, i)) * opts.binsize));
        end
        pn = nansum(pni);
        
        % Calculate P(x|n), log-probability of every position given observed spike vector
        pxn = pnx + px' - pn;
        
        % Convert to probability distribution over positions (i.e. sum to 1)
        pxn = exp(pxn) / nansum(exp(pxn));
        
    else

        pxn(1:size(fr_per_bin, 1), 1) = 0;

    end
    
end
