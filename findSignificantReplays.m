
function [sig_replay, loosely_sig_replay] = findSignificantReplays(gof, gof_cell_shuffled, gof_isi_shuffled, fr_per_bin, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts)

    % Find population bursts where goodness-of-fit exceeds all shuffles
    put_sig_replay = find(gof > cell2mat(cellfun(@(x, y) max([x y]), gof_cell_shuffled, gof_isi_shuffled, 'UniformOutput', false)));

    % Of this subset, find population bursts where absolute weighted correlation exceeds 0.5
    sig_replay = [];
    for i = 1:length(put_sig_replay)

        iBurst = put_sig_replay(i);

        % Find weighted correlation
        [weighted_corr(iBurst)] = weightedCorrelation(fr_per_bin, spiketimes, linearised_bin_centres, Timestamps_q, burst_start(iBurst), burst_end(iBurst), opts);

        % Include if absolute weighted correlation >= 0.5
        if weighted_corr(iBurst) >= 0.5 || weighted_corr(iBurst) <= -0.5
            sig_replay = [sig_replay; put_sig_replay(i)];
        end

end
    

    % Find loosely significant replay, where goodness-of-fit exceeds 95% of both null distributions
    put_loosely_sig_replay = find(gof > cell2mat(cellfun(@(x, y) max([prctile(x, 95) prctile(y, 95)]), gof_cell_shuffled, gof_isi_shuffled, 'UniformOutput', false)));
    put_loosely_sig_replay = put_loosely_sig_replay(~ismember(put_loosely_sig_replay, put_sig_replay));

    % Of this subset, find population bursts where absolute weighted correlation exceeds 0.5
    loosely_sig_replay = [];
    for i = 1:length(put_loosely_sig_replay)

        iBurst = put_loosely_sig_replay(i);

        % Find weighted correlation
        [weighted_corr(iBurst)] = weightedCorrelation(fr_per_bin, spiketimes, linearised_bin_centres, Timestamps_q, burst_start(iBurst), burst_end(iBurst), opts);

        % Include if absolute weighted correlation >= 0.5
        if weighted_corr(iBurst) >= 0.5 || weighted_corr(iBurst) <= -0.5
            loosely_sig_replay = [loosely_sig_replay; put_loosely_sig_replay(i)];
        end

    end

end