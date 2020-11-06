
function [] = plotCumulativeReplays(decoded_replay, decoded_replay_sh, pvalues, burst_start)

    figure; hold on

    % Plot shuffled replays
    mean_sh_replays = mean(cumsum(decoded_replay_sh, 2));
    std_sh_replays = std(cumsum(decoded_replay_sh, 2));
    patch([burst_start flip(burst_start)],...
        [mean_sh_replays + std_sh_replays flip(mean_sh_replays - std_sh_replays)],...
        'k', 'EdgeColor', 'none', 'FaceAlpha', 0.4)
    plot(burst_start, mean_sh_replays, 'k')

    % Plot actual replays
    plot(burst_start, cumsum(decoded_replay), 'k')
    hold on; scatter(burst_start(find(decoded_replay)), cumsum(decoded_replay(find(decoded_replay))), 'ko', 'filled')
    xlabel('Quiescence period')
    ylabel('Cumulative number of replays')

    % Plot p-values
    yyaxis right
    for i = 1:length(pvalues)-1
        if pvalues(i) < 0.05
            line([burst_start(i) burst_start(i+1)], [pvalues(i) pvalues(i+1)], 'Color', 'r')
        else
            line([burst_start(i) burst_start(i+1)], [pvalues(i) pvalues(i+1)], 'Color', 'b')
        end
    end
    i = pvalues < 0.05;
    scatter(burst_start(i), pvalues(i), 'ro', 'filled')
    i = pvalues >= 0.05;
    scatter(burst_start(i), pvalues(i), 'bo', 'filled')
    ax = gca; ax.YAxis(2).Color = 'r';
    ylim([0 1])
    ylabel('Monte Carlo p-value')

end