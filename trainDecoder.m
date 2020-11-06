
function [fr_per_bin, adaptive_bin_info] = trainDecoder(spiketimes, linearised_pos, linearised_bin_centres, Timestamps, run_starts, run_ends, opts)

    % Opts to update later: adaptive binning, convolving

    % Initialise running total of spikes per linearised bin per cell and occupancy
    bin_spikes_train = zeros(length(linearised_bin_centres), length(spiketimes));
    occupancy = zeros(length(linearised_bin_centres), 1);
    
    if opts.adaptive_binning
        
        if ~isfield(opts,'alpha') opts.alpha = 10.^6; end
        
        % Get timetamps and linearised position for all runs (also round timestamps to nearest milisecond)
        run_time = []; run_pos = [];
        for iRun = 1:length(run_starts)
            r_t = Timestamps(run_starts(iRun) : run_ends(iRun));
            run_time = [run_time round(r_t, 3)];
            run_pos = [run_pos linearised_pos(run_starts(iRun) : run_ends(iRun))];
        end
        
         % Get spiketimes for this run (also round spiketimes to nearest milisecond)
         run_spikes = [];
         for iRun = 1:length(run_starts)
            spk = cellfun(@(x) x(x >= Timestamps(run_starts(iRun)) & x <= Timestamps(run_ends(iRun))), spiketimes, 'UniformOutput', false);
            run_spikes = [run_spikes cellfun(@(x) round(x, 3), spk, 'UniformOutput', false)];
         end
         
         % Get position for every spike
         for iRun = 1:length(run_starts)
            [~, loc] = cellfun(@(x) ismember(x, run_time), run_spikes(:, iRun), 'UniformOutput', false);
            run_spikes_pos(iRun, :) = cellfun(@(x) run_pos(x), loc, 'UniformOutput', false);
         end
        
        adaptive_bin_info = []; % Record information about radius size etc.

        for iBin = 1:length(linearised_bin_centres)

            for iUnit = 1:length(spiketimes)

                r = 0; % Radius of adaptive bin
                centre = linearised_bin_centres(iBin); % Centre of adaptive bin
                done = false;

                while ~done

                    % Define circle around bin centre (ensure doesn't include positions beyond maze)
                    circle_bounds = [centre-r centre+r];
                    circle_bounds(circle_bounds < min(linearised_pos)) = min(linearised_pos);
                    circle_bounds(circle_bounds > max(linearised_pos)) = max(linearised_pos);
%                     circle_bounds = circle_bounds(circle_bounds >= min(linearised_pos) & circle_bounds <= max(linearised_pos));

                    % Find number of occupancy samples that fall within the circle
%                     N_occ = sum(ismember(run_pos, circle_bounds));
                    N_occ = sum(run_pos >= circle_bounds(1) & run_pos <= circle_bounds(2));

                    % Find number of spikes within the circle
%                     N_spikes = sum(ismember(cell2mat(run_spikes_pos(:, iUnit)'), circle_bounds));
                    N_spikes = sum(cell2mat(run_spikes_pos(:, iUnit)') >= circle_bounds(1) & cell2mat(run_spikes_pos(:, iUnit)') <= circle_bounds(2));

                    % Check if condition met
                    if (N_spikes*N_occ/1000) > opts.alpha / (N_occ.^2 * (r*100+0.5).^2) && N_occ > 0
                        fr_per_bin(iBin,iUnit) = (N_spikes/N_occ*1000);
                        done = true;
                    elseif range(circle_bounds) == range(linearised_pos) && N_occ > 0       % If true, have run out of maze and unit hasn't spiked enough
                        fr_per_bin(iBin,iUnit) = (N_spikes/N_occ*1000);
                        done = true;
                    elseif range(circle_bounds) == range(linearised_pos) && N_occ == 0     % If true, there is no position data and something is wrong
                        fr_per_bin(iBin,iUnit) = NaN;
                        done = true;
                    else
                        r = r + 0.01;
                    end

                end

                % Store info about adaptive binning
                adaptive_bin_info = [adaptive_bin_info; [iUnit iBin r N_occ N_spikes]];
            end

        end
        
    
    % Do it the easy way: spatial bins are uniform size
    else

        for iRun = 1:length(run_starts)

            % Get timetamps and linearised position for this run (also round timestamps to nearest milisecond)
            run_time = Timestamps(run_starts(iRun) : run_ends(iRun));
            run_time = round(run_time, 3);
            run_pos = linearised_pos(run_starts(iRun) : run_ends(iRun));

            % Get spiketimes for this run (also round spiketimes to nearest milisecond)
            run_spikes = cellfun(@(x) x(x >= Timestamps(run_starts(iRun)) & x <= Timestamps(run_ends(iRun))), spiketimes, 'UniformOutput', false);
            run_spikes = cellfun(@(x) round(x, 3), run_spikes, 'UniformOutput', false);

            % Get position for every spike
            [~, loc] = cellfun(@(x) ismember(x, run_time), run_spikes, 'UniformOutput', false);
            run_spikes_pos = cellfun(@(x) run_pos(x), loc, 'UniformOutput', false);

            % Get firing for every bin
            for iBin = 1:length(linearised_bin_centres)

                % Find number of spikes fired by each cell in this bin and add to running total
                bin_spikes_train(iBin,:) = bin_spikes_train(iBin,:) + cellfun(@(x) sum(x==linearised_bin_centres(iBin)), run_spikes_pos)';

                % Add occupancy to running total
                occupancy(iBin) = occupancy(iBin) + sum(run_pos == linearised_bin_centres(iBin));

            end

        end
        
        % Divide spikes by occupancy to get firing rate in Hz
        fr_per_bin = bin_spikes_train ./ repmat(occupancy, 1, length(spiketimes)) * 1000;
        
        adaptive_bin_info = NaN;
    
    end

    % Smooth
    try
        fr_per_bin = smoothts(fr_per_bin,'g',opts.smooth_kernel);
    catch
        fr_per_bin = smoothdata(fr_per_bin,'gaussian',opts.smooth_kernel);
    end

    % Plot
    if exist('opts.plot') && opts.plot==True
        figure; imagesc(fr_per_bin'); c = colorbar;
        xlabel('Linearised bin')
        ylabel('Pyramidal cells')
        c.Label.String = 'Mean firing rate (Hz)';
        drawnow
    end
            
end