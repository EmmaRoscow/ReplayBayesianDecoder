
function [decoded_duration, decoded_arm_coverage, decoded_weighted_corr, decoded_replay] = decodeWeightedCorr(fr_per_bin, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts)

    tic

    f = waitbar(0, opts.message);

    for iBurst = 1:length(burst_start)

        [decoded_duration(iBurst), decoded_arm_coverage(iBurst), decoded_weighted_corr(iBurst), decoded_replay(iBurst), pxn, traj_bins, traj_start, traj_end, bin_starts]...
            = decodeBurstWeightedCorr(fr_per_bin, spiketimes, linearised_bin_centres, Timestamps_q, burst_start(iBurst), burst_end(iBurst), opts);

        % Plot
        switch opts.plot
            
            case 'as_you_go'
                figure

                thresh = 1/size(pxn, 2) * 5;

                subplot(5, 1, 1:2)
                imagesc((bin_starts - bin_starts(1))*1000, linearised_bin_centres, pxn'); colormap(gca, hot); x = xlim;
                ylabel({'Linearised', 'bin position (m)'})
                c = colorbar; c.Label.String = 'Posterior prob.';

                subplot(5, 1, 3:4)
                imagesc((bin_starts - bin_starts(1))*1000, linearised_bin_centres, pxn' > thresh); colormap(gca, hot); x = xlim;
                ylabel({'Linearised', 'bin position (m)'})
                c = colorbar; c.Label.String = 'Posterior above threshold';
                c.Ticks = [0 1];
                c.TickLabels = {'No', 'Yes'};

                subplot(5, 1, 5); hold on
                line(x, [thresh thresh], 'Color', 'b')
                plot((bin_starts - bin_starts(1))*1000, max(pxn, [], 2), 'Color', [0.5 0.5 0.5]); xlim(x)
                plot((bin_starts(traj_bins(traj_start):traj_bins(traj_end)) - bin_starts(1)) * 1000, max(pxn(traj_bins(traj_start):traj_bins(traj_end), :), [], 2), 'Color', 'k', 'LineWidth', 1.5)
                xlabel('Time (ms)'); ylabel('MAP')

                drawnow
                
            case 'sig_only'
                
                if decoded_replay(iBurst)
                    
                    figure

                    thresh = 1/size(pxn, 2) * 5;

                    subplot(5, 1, 1:2)
                    imagesc((bin_starts - bin_starts(1))*1000, linearised_bin_centres, pxn'); colormap(gca, hot); x = xlim;
                    ylabel({'Linearised', 'bin position (m)'})
                    c = colorbar; c.Label.String = 'Posterior prob.';

                    subplot(5, 1, 3:4)
                    imagesc((bin_starts - bin_starts(1))*1000, linearised_bin_centres, pxn' > thresh); colormap(gca, hot); x = xlim;
                    ylabel({'Linearised', 'bin position (m)'})
                    c = colorbar; c.Label.String = 'Posterior above threshold';
                    c.Ticks = [0 1];
                    c.TickLabels = {'No', 'Yes'};

                    subplot(5, 1, 5); hold on
                    line(x, [thresh thresh], 'Color', 'b')
                    plot((bin_starts - bin_starts(1))*1000, max(pxn, [], 2), 'Color', [0.5 0.5 0.5]); xlim(x)
                    plot((bin_starts(traj_bins(traj_start):traj_bins(traj_end)) - bin_starts(1)) * 1000, max(pxn(traj_bins(traj_start):traj_bins(traj_end), :), [], 2), 'Color', 'k', 'LineWidth', 1.5)
                    xlabel('Time (ms)'); ylabel('MAP')

                    drawnow
                
                end
               
        end

        waitbar( iBurst / length(burst_start), f, {opts.message, [num2str(round(toc/60)) ' min(s) elapsed']})

    end

end