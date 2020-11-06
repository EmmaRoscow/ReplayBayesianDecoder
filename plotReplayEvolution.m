
function [] = plotReplayEvolution(Timestamps, t_start, t_end, Timestamps_q, burst_start, best_v, sig_replay, loosely_sig_replay, direction)


    % best_v > 0 means forward for rightward runs,
    % best_v < 0 means forward for leftward runs

    figure; hold on

    % Plot when the runs are
    for i = 1:length(t_start)
        patch(([Timestamps(t_start(i)) Timestamps(t_end(i)) Timestamps(t_end(i)) Timestamps(t_start(i))] - sessInfo.Epochs.MazeEpoch(1))/60, [0 0 1 1], 'k')
    end

    % Plot population bursts
    scatter((Timestamps_q(burst_start) - sessInfo.Epochs.MazeEpoch(1))/60, ones(1, length(burst_start))*1.5, 5, 'g', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.25)

    % Plot forward and backward significant replays
    switch direction
        case "rightward"
            forward_replay = find(best_v(sig_replay) > 0);
            backward_replay = find(best_v(sig_replay) < 0);
        case "leftward"
            forward_replay = find(best_v(sig_replay) < 0);
            backward_replay = find(best_v(sig_replay) > 0);
    end
    scatter((Timestamps_q(burst_start(sig_replay(forward_replay))) - sessInfo.Epochs.MazeEpoch(1))/60, ones(1, length(forward_replay))*2, 'b', 'filled', 'MarkerFaceAlpha', 0.25)
    scatter((Timestamps_q(burst_start(sig_replay(backward_replay))) - sessInfo.Epochs.MazeEpoch(1))/60, ones(1, length(backward_replay))*3, 'm', 'filled', 'MarkerFaceAlpha', 0.25)

    % Plot forward and backward loosely significant replays
    switch direction
        case "rightward"
            forward_replay = find(best_v(loosely_sig_replay) > 0);
            backward_replay = find(best_v(loosely_sig_replay) < 0);
        case "leftward"
            forward_replay = find(best_v(sig_replay) < 0);
            backward_replay = find(best_v(sig_replay) > 0);
    end
    scatter((Timestamps_q(burst_start(loosely_sig_replay(forward_replay))) - sessInfo.Epochs.MazeEpoch(1))/60, ones(1, length(forward_replay))*2.25, 'b', 'filled', 'MarkerFaceAlpha', 0.25)
    scatter((Timestamps_q(burst_start(loosely_sig_replay(backward_replay))) - sessInfo.Epochs.MazeEpoch(1))/60, ones(1, length(backward_replay))*3.25, 'm', 'filled', 'MarkerFaceAlpha', 0.25)

    set(gca, 'YTick', [0.5 1.5 2 2.25 3 3.25])
    set(gca, 'YTickLabel', {'Runs', 'Pop. bursts', 'Forward replays (sig)', 'Forward replays (loosely sig)', 'Backward replays (sig)', 'Backward replays (loosely sig)'})
    xlabel('Time (mins)')

end