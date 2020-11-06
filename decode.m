
function [pxn, max_gof, best_v, best_c, n, bin_starts, bin_ends] = decode(fr_per_bin, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts)

    % Divide popluation burst into bins of specified length, with specified step size
    bin_starts = round(Timestamps_q(burst_start),3) : opts.stepsize : round(Timestamps_q(burst_end),3)-(opts.binsize-1)/1000;
    bin_ends = bin_starts + opts.binsize;
    
    %% Initialise for decoding
    n = []; pn = []; pnx = []; pxn = [];
    
    % Get number of spikes
%     bin_starts = round(bin_starts, 3);
%     bin_ends = round(bin_ends, 3);
%     n_test = cellfun(@(spk, x, y) sum(round(spk, 3) >= x & round(spk, 3) < y), repmat(spiketimes, 1, length(bin_starts)), repmat(mat2cell(bin_starts, 1, ones(81, 1)), length(spiketimes), 1), repmat(mat2cell(bin_ends, 1, ones(81, 1)), length(spiketimes), 1))';
    
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
    
    
    % Fit lines
    velocities = [-3 -2 -1 1 2 3]; % If spatial bin is 3cm and stepsize is 2ms, 3 bins/step is 45 metres/second
    start_pos = -1 * size(pxn, 2) : 2 * size(pxn, 2);
    gof = fitLines(pxn, velocities, start_pos, bin_starts, linearised_bin_centres, opts);
    
    % Find the best fit
    max_gof = max(gof(:));
    [best_c, best_v] = find(gof==max_gof, 1);
    best_v = velocities(best_v);
    best_c = start_pos(best_c);
    
    % Plot
    switch opts.plot
        
        case {'as_you_go', 'best_only'}
        
            % Calculate again
            y  = best_v*([1:size(pxn, 1)]-1) + best_c;
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
            subplot(3, 1, 1:2); hold on
            imagesc(bin_starts - bin_starts(1), linearised_bin_centres, pxn', [0 1])
            colormap(gca, hot)
            set(gca, 'YDir', 'normal')
            set(gca, 'XTickLabel', [])
            ylabel('Linearised position (m)')
            title(['Line of best fit: ' num2str(max_gof)])
            axis tight; x = xlim;
            line(bin_starts(pos) - bin_starts(1), linearised_bin_centres(cellfun(@(x) min(x), band_idx)), 'Color', 'w')
            line(bin_starts(pos) - bin_starts(1), linearised_bin_centres(cellfun(@(x) max(x), band_idx)), 'Color', 'w')
            scatter(bin_starts(pos) - bin_starts(1), linearised_bin_centres(y), 8, 'ow', 'filled')
            subplot(3, 1, 3)
            plot(bin_starts(pos) - bin_starts(1), gof_per_bin)
            xlim(x)
            ylim([0 1])
            xlabel('Time bin')
            ylabel('Goodness-of-fit')
            drawnow
            
    end
        
end