
function [] = plotSignificantReplay(to_plot, fr_per_bin, best_v, best_c, gof, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts)

    for iReplay  = 1:length(to_plot)

        iBurst = to_plot(iReplay);

        % Divide popluation burst into bins of specified length, with specified step size
        bin_starts = round(Timestamps_q(burst_start(iBurst)),3) : opts.stepsize : round(Timestamps_q(burst_end(iBurst)),3)-(opts.binsize-1)/1000;
        bin_ends = bin_starts + opts.binsize;

        % Initialise for decoding
        n = []; pn = []; pnx = []; pxn = [];

        for t = 1:length(bin_starts)

            % Get number of spikes
            time_bin_start = round(bin_starts(t),3);
            time_bin_end = round(bin_ends(t),3);
            n(t,:) = cellfun(@(x) sum(x >= bin_starts(t) & x < bin_ends(t)), spiketimes);

            % Check number of active cells is at least 4
            if sum(n(t,:) > 0) >= 4

                % Calculate P(n|x), log-probability of spikes assuming Poisson distribution for each binned position
                for i = 1:length(spiketimes)
                    pnix(:,i) = log(pdf('Poisson', n(t,i), fr_per_bin(:,i) * opts.binsize));
                end
                pnix(isinf(pnix)) = -100; % Cap low probabilities
                pnx(t,:) = sum(pnix,2);

                % Calculate P(x), log-probability of every position; assume uniform
                px = log(ones(1,size(fr_per_bin,1))/size(fr_per_bin,1));

                % Calculate P(n), log-probability of observed spike vector; assume Poisson firing rates
                for i = 1:length(spiketimes)
                    pni(i) = log(pdf('Poisson', n(t, i), nanmean(fr_per_bin(:, i)) * opts.binsize));
                end
    %             pni(isinf(pni)) = -100; % Cap low probabilities
                pn(t) = nansum(pni);

                % Calculate P(x|n), log-probability of every position given observed spike vector
                pxn(t,:) = pnx(t,:)'+px'-pn(t);

                % Convert to probability distribution over positions (i.e. sum to 1)
                pxn(t,:) = exp(pxn(t, :)) / nansum(exp(pxn(t, :)));

            else

                pxn(t,1:length(linearised_bin_centres)) = 0;

            end

        end

        % Calculate again
        y  = best_v(iBurst)*([1:size(pxn, 1)]-1) + best_c(iBurst);
        pos = find(y >= 1 & y <= size(pxn, 2));
        y = y(pos);
        gof_per_bin = zeros(1, length(y));
        band_idx = {}; gof_per_bin = [];
        for t = 1:length(y)
            band_idx{t} = y(t)-4 : y(t)+4;
            band_idx{t} = band_idx{t}(band_idx{t} >=1 & band_idx{t} <= size(pxn, 2));
            gof_per_bin(t) = nansum(pxn(pos(t), band_idx{t})) / nansum(pxn(pos(t), :));
        end
        gof_per_bin(isnan(gof_per_bin)) = 0;

         % Plot
        figure
        subplot(10, 2, 1:8); hold on
        imagesc(bin_starts - bin_starts(1), linearised_bin_centres, pxn')
        axis tight
        colormap(gca, hot)
        set(gca, 'YDir', 'normal')
        set(gca, 'XTickLabel', [])
        ylabel({'Linearised', 'position (m)'})
        title(['Line of best fit: ' num2str(round(gof(iBurst), 2)) ' // weighted correlation: ' num2str(round(weighted_corr(iBurst), 2))])
        axis tight; x = xlim;
        line(bin_starts(pos) - bin_starts(1), linearised_bin_centres(cellfun(@(x) min(x), band_idx)), 'Color', 'w')
        line(bin_starts(pos) - bin_starts(1), linearised_bin_centres(cellfun(@(x) max(x), band_idx)), 'Color', 'w')
        % scatter(bin_starts(pos) - bin_starts(1), linearised_bin_centres(y), 8, 'ow', 'filled')

        subplot(10, 2, 9:12)
        plot(bin_starts(pos) - bin_starts(1), gof_per_bin, 'k')
        xlim(x)
        ylim([0 1])
        xlabel('Time bin start (s)')
        ylabel({'Goodness-', 'of-fit'})

        subplot(10, 2, [17 19])
        histogram(gof_cell_shuffled{iBurst}(:), 0:0.01:1, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none')
        y = ylim; xlim([0 1])
        line([gof(iBurst) gof(iBurst)], y, 'Color', 'b')
        xlabel('Goodness of fit');
        ylabel({'Null', 'distribution'})
        title('Shuffling cell identities')

        subplot(10, 2, [18 20])
        histogram(gof_isi_shuffled{iBurst}(:), 0:0.01:1, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none')
        y = ylim; xlim([0 1])
        line([gof(iBurst) gof(iBurst)], y, 'Color', 'b')
        xlabel('Goodness of fit');
        title('Shuffling interspike intervals')    
        l = legend('Null dist''n', 'Actual g.o.f.'); set(l, 'Box', 'off');

        drawnow

    end
    
end