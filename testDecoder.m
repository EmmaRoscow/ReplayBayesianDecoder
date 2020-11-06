
function [position_error, n, est_binned_pos, true_binned_pos, pxn] = testDecoder(fr_per_bin, spiketimes, linearised_pos, linearised_bin_centres, Timestamps, run_start, run_end, opts)

    % Divide test run into bins of specified length, with specified step size
    bin_starts = round(Timestamps(run_start),3) : opts.stepsize : round(Timestamps(run_end),3)-(opts.binsize-1)/1000;
    bin_ends = bin_starts + opts.binsize;
        
    % Initialise for decoding
    est_binned_pos = [];
    true_binned_pos = [];
    n = []; pn = []; pnx = []; pxn = [];
    counter = 1;
    
    f = waitbar(0, 'Beginning decoding...'); tic
    
    for t = 1:length(bin_starts)

        % Get number of spikes
        time_bin_start = round(bin_starts(t),3);
        time_bin_end = round(bin_ends(t),3);
        n(counter,:) = cellfun(@(x) sum(round(x,3) >= bin_starts(t) & round(x,3) < bin_ends(t)), spiketimes);

        % Check number of active cells is at least 4
        if sum(n(counter,:) > 0) >= 4

            % Calculate P(n|x), log-probability of spikes assuming Poisson distribution for each binned position
            for i = 1:length(spiketimes)
                pnix(:,i) = log(pdf('Poisson',n(counter,i),fr_per_bin(:,i) * opts.binsize));
            end
            pnix(isinf(pnix)) = -100; % Cap low probabilities
            pnx(counter,:) = sum(pnix,2);
%             pnx(isnan(pnx)) = -inf;

            % Calculate P(x), log-probability of every position; assume
            % uniform
            px = log(ones(1,size(fr_per_bin,1))/size(fr_per_bin,1));
            
            % Calculate P(n), log-probability of observed spike vector
            px = log(ones(1,size(fr_per_bin,1))/size(fr_per_bin,1));
            pn(counter) = sum(pnx(counter,:)+px);
            
            % Calculate P(n), log-probability of observed spike vector;
            % assume Poisson firing rates
            for i = 1:length(spiketimes)
                pni(i) = log(pdf('Poisson', n(counter, i), nanmean(fr_per_bin(:, i)) * opts.binsize));
            end
            pn(counter) = nansum(pni);

            % Calculate P(x|n), log-probability of every position given
            % observed spike vector
%             nBins = size(fr_per_bin,1);
%             pxn(counter,1:nBins) = smooth(pnx(counter,1:nBins)'+px(1:nBins)'-pn(counter),5);
%             pxn(counter,:) = pnx(counter,:)'+px'-pn(counter);
            pxn(counter, :) = pnx(counter, :)' + px' - pn(counter);
            if ~isfield(opts, "normalised") || opts.normalised == True
                pxn(counter, :) = exp(pxn(counter, :)) / nansum(exp(pxn(counter, :)));
            end

            % Decode position
%             [prob,estimated_bin] = max(pnx(counter,:));
            [prob, estimated_bin] = max(pxn(counter, :));
            est_binned_pos(counter) = linearised_bin_centres(estimated_bin);

            % Calculate error as distance from true position
            idx1 = find(round(Timestamps,3) == round(time_bin_start,3));
            idx2 = find(round(Timestamps,3) == round(time_bin_end,3));
            true_binned_pos(counter) = mean(linearised_pos(idx1 : idx2));
            position_error(counter) = abs(true_binned_pos(counter) - est_binned_pos(counter));

            % Update counter
            counter = counter + 1;

        end
        
        waitbar(t / length(bin_starts), f, {['Decoded ' num2str(t) ' of ' num2str(length(bin_starts)) ' time bins, '], [num2str(round(toc, 1)) ' seconds']})
        
    end
    
    % Plot
    if isfield(opts, 'plot') && opts.plot == true
        figure('Position', [680 483 560 495])
        subplot(3, 2, [4 6]); hold on
        plot(true_binned_pos, 'k', 'LineWidth', 2)
        plot(est_binned_pos, 'b', 'LineWidth', 2)
        xlabel('Time bin')
        ylabel('Linearised position')
        l = legend('True position', 'Decoded (estimated) postiion', 'Location', 'southoutside');
        l.Box = 'off';
        xlim([1 length(est_binned_pos)])
        y = ylim;
        
        subplot(3, 2, [1 3])
        imagesc(fr_per_bin'); colorbar; colormap parula
        ylabel('Pyramidal cell #')
        title({'Linearised firing rate maps', '(training data)'})
        
        subplot(3, 2, 2)
%         norm_pnx = exp(pnx) ./ repmat(sum(exp(pnx'))', 1, size(pnx, 2));
%         imagesc(1:size(pxn, 1), linearised_bin_centres, norm_pnx)
        imagesc(1:size(pnx, 1), linearised_bin_centres, pxn'); colormap(gca, hot)
        axis tight
        set(gca, 'YDir', 'normal')
%         xlabel('Time bin')
        set(gca, 'XTickLabel', [])
        ylabel('Linearised position')
        title('Estimated probability')
        xlim([1 length(est_binned_pos)])
        
        subplot(3, 2, 5)
        imagesc(n'); colormap parula
        colorbar
        xlabel('Linearised position')
        title('Observed firing rate (test data)')
    end
        
    close(f)

end