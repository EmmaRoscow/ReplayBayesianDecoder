
function [gof] = decodeSignificance(fr_per_bin, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, best_c, best_v, opts)

    % Divide population burst into bins of specified length, with specified step size
    bin_starts = round(Timestamps_q(burst_start),3) : opts.stepsize : round(Timestamps_q(burst_end),3)-(opts.binsize-1)/1000;
    bin_ends = bin_starts + opts.binsize;
    
    % Shifting involves changing spiketimes and/or fr_per_bin, so record the original values here
    spiketimes_ = spiketimes;
    fr_per_bin_ = fr_per_bin;
    
    for iSh = 1:opts.nSh
    
        switch opts.shuffle
            
            case 'cell_identities'
            % Shuffle cell identities
            reorder = randperm(length(spiketimes_));
            spiketimes = spiketimes_(reorder);

            case 'interspike_intervals'
            % Shuffle interspike intervals
            spiketimes_ = cellfun(@(x) x(x >= bin_starts(1) & x < bin_ends(end)), spiketimes_, 'UniformOutput', false);
            isi = cellfun(@(x) diff([bin_starts(1); x]), spiketimes_, 'UniformOutput', false);
            isi = cellfun(@(x) x(randperm(length(x))), isi, 'UniformOutput', false);
            spiketimes = cellfun(@(x) bin_starts(1) + cumsum(x)', isi, 'UniformOutput', false);
            
            case 'circular_shift'
            % Circularly shift each cell's firing rate map
            shift = randi(size(fr_per_bin_, 1), size(fr_per_bin_, 2), 1);
            reorder = rem([1:size(fr_per_bin_, 1)] + shift - 1, size(fr_per_bin_, 1)) + 1;
            fr_per_bin = fr_per_bin_(reorder');
            
        end
   
        % Do decoding
        pxn = arrayfun(@(x, y) decodeBin(x, y, spiketimes, fr_per_bin, opts), bin_starts, bin_ends, 'UniformOutput', false);
        pxn = cell2mat(pxn)';

%         for t = 1:length(bin_starts)
%             pxn(t, :) = decodeBin(round(bin_starts(t), 3), round(bin_ends(t), 3), spiketimes, fr_per_bin, opts);
%         end
        
        % Fit lines
        velocities = [-3 -2 -1 1 2 3]; % If spatial bin is 3cm and stepsize is 2ms, 3 bins/step is 45 metres/second
        start_pos = -1 * size(pxn, 2) : 2 * size(pxn, 2);
        gof_iSh = fitLines(pxn, velocities, start_pos, bin_starts, linearised_bin_centres, opts);
        gof(iSh) = max(gof_iSh(:));
                
    end
   
end