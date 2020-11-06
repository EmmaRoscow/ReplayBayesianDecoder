
function [] = plotJointReplays(fr_per_bin_leftward, fr_per_bin_rightward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts)

    for iBurst = 1:length(burst_start)

        figure

        % Decode rightward
        [decoded_duration, decoded_arm_coverage, decoded_weighted_corr, decoded_replay, pxn_right, traj_bins, traj_start, traj_end, bin_starts]...
                = decodeBurstWeightedCorr(fr_per_bin_rightward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start(iBurst), burst_end(iBurst), opts);

        thresh = 1/size(pxn_right, 2) * 5;

        subplot(6, 1, 1:2)
        imagesc((bin_starts - bin_starts(1))*1000, linearised_bin_centres, pxn_right'); colormap(gca, hot); x = xlim;
        ylabel({'Linearised', 'bin position (m)'})
        c = colorbar; c.Label.String = 'Posterior prob.';
        title('Rightward decoding')
        pos = get(gca, 'Position');
        x = xlim;
        set(gca, 'XTick', [])

    %     subplot(6, 1, 3)
    %     line(x, [thresh thresh], 'Color', 'b')
    %     plot((bin_starts - bin_starts(1))*1000, max(pxn, [], 2), 'Color', [0.5 0.5 0.5]); xlim(x)
    %     plot((bin_starts(traj_bins(traj_start):traj_bins(traj_end)) - bin_starts(1)) * 1000, max(pxn(traj_bins(traj_start):traj_bins(traj_end), :), [], 2), 'Color', 'k', 'LineWidth', 1.5)
    %     xlabel('Time (ms)'); ylabel('MAP')
    %     pos2 = get(gca, 'Position'); set(gca, 'Position', [pos2(1) pos2(2) pos(3) pos2(4)]);
    %     xlim(x);

        % Decode leftward
        [decoded_duration, decoded_arm_coverage, decoded_weighted_corr, decoded_replay, pxn_left, traj_bins, traj_start, traj_end, bin_starts]...
                = decodeBurstWeightedCorr(fr_per_bin_leftward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start(iBurst), burst_end(iBurst), opts);

        thresh = 1/size(pxn_left, 2) * 5;

        subplot(6, 1, 3:4)
        imagesc((bin_starts - bin_starts(1))*1000, linearised_bin_centres, pxn_left'); colormap(gca, hot); x = xlim;
        ylabel({'Linearised', 'bin position (m)'})
        c = colorbar; c.Label.String = 'Posterior prob.';
        title('Leftward decoding')
        pos = get(gca, 'Position');
        x = xlim;

    %     subplot(6, 1, 6)
    %     line(x, [thresh thresh], 'Color', 'b')
    %     plot((bin_starts - bin_starts(1))*1000, max(pxn, [], 2), 'Color', [0.5 0.5 0.5]); xlim(x)
    %     plot((bin_starts(traj_bins(traj_start):traj_bins(traj_end)) - bin_starts(1)) * 1000, max(pxn(traj_bins(traj_start):traj_bins(traj_end), :), [], 2), 'Color', 'k', 'LineWidth', 1.5)
    %     xlabel('Time (ms)'); ylabel('MAP')
    %     pos2 = get(gca, 'Position'); set(gca, 'Position', [pos2(1) pos2(2) pos(3) pos2(4)]);
    %     xlim(x);

        subplot(6, 1, 5:6); hold on
        line(x, [0 0], 'Color', 'k', 'LineStyle', ':')
        cc = []; rl = []; ru = [];
        for i = 1:size(pxn_right, 1)
            [ans1, ~, ans2, ans3] = corrcoef(pxn_right(i, :), pxn_left(i, :));
            cc(i) = ans1(2);
            rl(i) = ans2(2);
            ru(i) = ans3(2);
        end
        cc(isnan(cc)) = 0; rl(isnan(rl)) = 0; ru(isnan(ru)) = 0;
        patch([(bin_starts - bin_starts(1))*1000 flip((bin_starts - bin_starts(1))*1000)], [rl flip(ru)], 'k', 'FaceAlpha', 0.25, 'EdgeColor', 'none')
        plot((bin_starts - bin_starts(1))*1000, cc, 'k')
        pos2 = get(gca, 'Position'); set(gca, 'Position', [pos2(1) pos2(2) pos(3) pos2(4)]);
        xlim(x); ylim([-1.1 1.1]); box off
        ylabel('Correlation')
        xlabel('Time (ms)')
        
        drawnow
    
    end
    
end