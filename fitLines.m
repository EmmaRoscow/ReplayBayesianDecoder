function [gof] = fitLines(pxn, velocities, start_pos, bin_starts, linearised_bin_centres, opts)
    
    % Plot
    switch opts.plot case 'as_you_go'
            figure; subplot(10, 4, [1:4:40 2:4:40]); hold on
            imagesc(bin_starts - bin_starts(1), linearised_bin_centres, pxn')
            axis tight
            set(gca, 'YDir', 'normal')
            ylabel('Linearised position (m)')
            xlabel('Time bin')
            title('Estimated probability, P(x|n)')
            iSub = 3;
            drawnow
    end
        
    for c = 1:length(start_pos)

        for v = 1:length(velocities)
            
            V = velocities(v);
            C = start_pos(c);
            
            y = V*([1:size(pxn, 1)]-1) + C;
            pos = find(y >= 1 & y <= size(pxn, 2));
            y = y(pos);
            
            % If this line doesn't fit the burst at all, skip
            if isempty(y)
                continue
            end

            gof_per_bin = zeros(1, length(y));
            band_idx = {}; gof_per_bin = [];

            % Identify which spatial bins fall into allowed band
            for t = 1:length(pos)
                band_idx{t} = y(t)-4 : y(t)+4;
                band_idx{t} = band_idx{t}(band_idx{t} >=1 & band_idx{t} <= size(pxn, 2));
            end
            
            % Plot
            switch opts.plot case 'as_you_go'
                    subplot(10, 4, [1:4:40 2:4:40]);
                    l1 = line(bin_starts(pos) - bin_starts(1), linearised_bin_centres(cellfun(@(x) min(x), band_idx)), 'Color', 'w');
                    l2 = line(bin_starts(pos) - bin_starts(1), linearised_bin_centres(cellfun(@(x) max(x), band_idx)), 'Color', 'w');
                    sc = scatter(bin_starts(pos) - bin_starts(1), linearised_bin_centres(y), 8, 'ow', 'filled');
                    drawnow
            end

            % Find proportion of probability distribution in each band
            for t = 1:length(pos)
                gof_per_bin(t) = nansum(pxn(pos(t), band_idx{t})) / nansum(pxn(pos(t), :));
            end

            % Find overall goodness of fit
            gof_per_bin(isnan(gof_per_bin)) = 0;
            gof(c, v) = sum(gof_per_bin) /  size(pxn, 1);
            
            % Plot
            switch opts.plot case 'as_you_go'
                delete(l1); delete(l2); delete(sc);
                if iSub <= 40
                    subplot(10, 4, iSub)
                    imagesc(bin_starts - bin_starts(1), linearised_bin_centres, pxn')
                    set(gca, 'YDir', 'normal')
                    line(bin_starts(pos) - bin_starts(1), linearised_bin_centres(cellfun(@(x) min(x), band_idx)), 'Color', 'w');
                    line(bin_starts(pos) - bin_starts(1), linearised_bin_centres(cellfun(@(x) max(x), band_idx)), 'Color', 'w');
                    title(num2str(gof(c, v)), 'FontSize', 8)
                    set(gca, 'XTick', [])
                    if iSub == 39 iSub = 4;
                    else iSub = iSub + 4; end
                end
                drawnow
            end
        
            
        end
        
    end
    
end