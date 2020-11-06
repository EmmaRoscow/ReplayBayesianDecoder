
function [decoded_duration, decoded_arm_coverage, decoded_weighted_corr, decoded_replay, decoded_replay_combined] = permutationWeightedCorr(fr_per_bin_1, fr_per_bin_2, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts)

    decoded_duration_1 = []; decoded_duration_2 = [];
    decoded_arm_coverage_1 = []; decoded_arm_coverage_2 = [];
    decoded_weighted_corr_1 = []; decoded_weighted_corr_2 = [];
    decoded_replay_1 = []; decoded_replay_2 = [];
    spiketimes__ = spiketimes;
    fr_per_bin_1_ = fr_per_bin_1;
    fr_per_bin_2_ = fr_per_bin_2;

    f = waitbar(0, opts.message);

    tic

    for iSh = 1:opts.nSh

        switch opts.shuffle
            case 'circular_shifts'
            % Circularly shift rate maps
            shift = randi(size(fr_per_bin_, 1), size(fr_per_bin_, 2), 1);
            reorder = rem([1:size(fr_per_bin_, 1)] + shift - 1, size(fr_per_bin_, 1)) + 1;
            fr_per_bin_1 = fr_per_bin_1_(reorder');
            fr_per_bin_2 = fr_per_bin_2_(reorder');
        end

        % Decode replay
        parfor iBurst = 1:length(burst_start)
            
            spiketimes = [];

            switch opts.shuffle
                case 'interspike_interval'
                % Shuffle interspike intervals
                spiketimes_ = cellfun(@(x) x(x >= Timestamps_q(burst_start(iBurst)) & x < Timestamps_q(burst_end(iBurst))), spiketimes__, 'UniformOutput', false);
                isi = cellfun(@(x) diff([Timestamps_q(burst_start(iBurst)); x]), spiketimes_, 'UniformOutput', false);
                isi = cellfun(@(x) x(randperm(length(x))), isi, 'UniformOutput', false);
                spiketimes = cellfun(@(x) Timestamps_q(burst_start(iBurst)) + cumsum(x)', isi, 'UniformOutput', false);
            end

            % Decode for first direction
            [decoded_duration_1(iSh, iBurst), decoded_arm_coverage_1(iSh, iBurst), decoded_weighted_corr_1(iSh, iBurst), decoded_replay_1(iSh, iBurst)]...
            = decodeBurstWeightedCorr(fr_per_bin_1, spiketimes, linearised_bin_centres, Timestamps_q, burst_start(iBurst), burst_end(iBurst), opts);
        
            % Decode for second direction
            [decoded_replay_2(iSh, iBurst), decoded_arm_coverage_2(iSh, iBurst), decoded_weighted_corr_2(iSh, iBurst), decoded_replay_2(iSh, iBurst)]...
            = decodeBurstWeightedCorr(fr_per_bin_2, spiketimes, linearised_bin_centres, Timestamps_q, burst_start(iBurst), burst_end(iBurst), opts);
        
        end

        waitbar(iSh/opts.nSh, f, {opts.message, [num2str(round(toc/60)) ' min(s) elapsed']})

    end
    
    % Combine
    decoded_duration = [decoded_duration_1; decoded_duration_2];
    decoded_arm_coverage = [decoded_arm_coverage_1; decoded_arm_coverage_2];
    decoded_weighted_corr = [decoded_weighted_corr_1; decoded_weighted_corr_2];
    decoded_replay = [decoded_replay_1; decoded_replay_2];
    decoded_replay_combined = decoded_replay_1 + decoded_replay_2;

end