        
function [decoded_duration, decoded_arm_coverage, decoded_weighted_corr, decoded_replay, pxn, traj_bins, traj_start, traj_end, bin_starts] = decodeBurstWeightedCorr(fr_per_bin, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts)

    % Get pxn again
    [pxn, ~, ~, ~, ~, bin_starts, bin_ends] = decode(fr_per_bin, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts);

    % Smooth
    pxn = smoothdata(pxn, 'gaussian', 5);

    % See which bins exceed threshold
    thresh = 1/size(pxn, 2) * 5;
    traj_bins = find(max(pxn, [], 2) > thresh);

    % Check decoded decoded_traj meets other criteria
    if isempty(traj_bins)
        decoded_duration = 0;
        decoded_arm_coverage = 0;
        decoded_weighted_corr = 0;
        decoded_replay = 0;
        traj_start = [];
        traj_end = [];

    else
        
        % Find largest peak of maximum a priori
        map = max(pxn(traj_bins, :), [], 2);
        [~, peak] = max(map);

        % Find beginning and end of above-threshold trajectory around peak
        traj_start = find(diff(traj_bins(1:peak)) > 1, 1, 'last') + 1;
        traj_end = find(diff(traj_bins(peak:end)) > 1, 1, 'first') + peak - 1;
        if isempty(traj_start) traj_start = 1; end
        if isempty(traj_end) traj_end = length(traj_bins); end

        % Check duration
        decoded_duration = bin_ends(traj_bins(traj_end)) - bin_starts(traj_bins(traj_start));

        % Check arm coverage
        decoded_arm_coverage = sum(max(pxn(traj_bins(traj_start):traj_bins(traj_end), :), [], 1) > thresh) / size(pxn, 2);

        % Check absolute weighted correlation
        decoded_weighted_corr = weightedCorrelation(pxn(traj_bins(traj_start):traj_bins(traj_end), :));

        % Check if meets criteria for inclusion as a replay
        if decoded_duration >= 0.05 && decoded_arm_coverage >= 0.1 && abs(decoded_weighted_corr) > 0.5
            decoded_replay = true;
        else
            decoded_replay = false;
        end
        
    end
    
end