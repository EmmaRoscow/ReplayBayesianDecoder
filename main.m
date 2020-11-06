
% Filepaths and parameters to change
filepath = 'C:\Users\Emma Roscow\Documents\GrosmarkBuzsàki data\NoveltySessInfoMatFiles';
filename = 'Cicero_09012014_sessInfo.mat';
opts.save_name = [filename(1:find(filename == '.')-1) '_data'];

% Load data for this session
load([filepath filesep filename])
addpath('C:\Users\Emma Roscow\Documents\GrosmarkBuzsàki data\Bayesian decoding')

%% ===== Set up ===== %%

% Define environment for binning
binsize = 0.03; % In metres

% Restrict position data to maze epoch
Xpos = sessInfo.Position.TwoDLocation(:, 1);
Ypos = sessInfo.Position.TwoDLocation(:, 2);
Timestamps = sessInfo.Position.TimeStamps;
[start_index, end_index] = truncateData(filename, Timestamps, Xpos, Ypos);

[Xpos, Ypos, Timestamps] = processPositionData(Xpos, Ypos, Timestamps, sessInfo.Epochs.MazeEpoch, start_index, end_index);

% Linearise positions
[linearised_pos, linearised_bin_centres, pos_bins] = linearisePositions(Xpos, Ypos, binsize, sessInfo.Position.MazeType);

% Get speed during maze epoch
kernel = 1000;
speed = getSpeed(Xpos, Ypos, Timestamps, kernel);

% Get rightward runs
speed_thresh = 0.05;
[leftward_t_start, leftward_t_end, rightward_t_start, rightward_t_end] = getRuns(linearised_pos, speed, speed_thresh, Timestamps, Xpos, opts);

% Get spiketimes for pyramidal cells
spiketimes = getSpiketimes(sessInfo.Spikes.SpikeTimes, sessInfo.Spikes.SpikeIDs, sessInfo.Spikes.PyrIDs);

%% ===== Find population bursts during maze epoch ===== %%

[burst_start, burst_end, quiescence_mua, Timestamps_q] = getPopulationBursts(Timestamps, sessInfo.Spikes.SpikeTimes, speed);

%% ===== Train and test decoder ===== %%

% Test decoder (if desired): train on every run but one and compare predicted position to actual position
opts.binsize = 0.15;
opts.stepsize = 0.05;
opts.smooth_kernel = 3;
opts.plot = true;
opts.adaptive_binning = true;
opts.spatial_bin_size = binsize;
for iRun = 1:length(rightward_t_start)
    
    % Train on all other runs
    opts.binsize = 0.05; opts.stepsize = 0.005;
    fr_per_bin_rightward = trainDecoder(spiketimes, linearised_pos, linearised_bin_centres, Timestamps, rightward_t_start([1:iRun-1 iRun+1:end]), rightward_t_end([1:iRun-1 iRun+1:end]), opts);
    
    % Test on odd-run-out
    opts.binsize = 0.15; opts.stepsize = 0.05;
    [position_error_rightward{iRun}, n{iRun}, est_binned_pos_rightward{iRun}, true_binned_pos_rightward{iRun}] = testDecoder(fr_per_bin_rightward, spiketimes, linearised_pos, linearised_bin_centres, Timestamps, rightward_t_start(iRun), rightward_t_end(iRun), opts);
    
    set(gcf, 'Name', num2str(iRun))
    
    fprintf('%s Completed test for rightward run %i, median position error %.1f cm\n', datetime('now'), iRun, median(position_error_rightward{iRun})*100)
    
end

for iRun = 1:length(leftward_t_start)
    
    % Train on all other runs
    opts.binsize = 0.05; opts.stepsize = 0.005;
    fr_per_bin_leftward = trainDecoder(spiketimes, linearised_pos, linearised_bin_centres, Timestamps, leftward_t_start([1:iRun-1 iRun+1:end]), leftward_t_end([1:iRun-1 iRun+1:end]), opts);
    
    % Test on odd-run-out
    opts.binsize = 0.15; opts.stepsize = 0.05;
    [position_error_leftward{iRun}, n{iRun}, est_binned_pos_leftward{iRun}, true_binned_pos_leftward{iRun}] = testDecoder(fr_per_bin_leftward, spiketimes, linearised_pos, linearised_bin_centres, Timestamps, leftward_t_start(iRun), leftward_t_end(iRun), opts);
    
    set(gcf, 'Name', num2str(iRun))
    
    fprintf('%s Completed test for leftward run %i, median position error %.1f cm\n', datetime('now'), iRun, median(position_error_leftward{iRun})*100)
    
end
    
% Plot position errors
figure; hold on
for i = 1:length(position_error_rightward)
    xjitter = linspace(i-0.4, i+0.4, length(position_error_rightward{i}));
    scatter(xjitter(randperm(length(xjitter))), position_error_rightward{i}*100, '.', 'MarkerEdgeColor', [0.25 0.25 0.25], 'MarkerEdgeAlpha', 0.25)
%     boxplot(position_error{i}, 'Colors', 'k', 'Widths', 0.5, 'Whisker', inf, 'Position', i);
end
plot(cellfun(@(x) median(x)*100, position_error_rightward), 'b', 'LineWidth', 2)
xlim([0 length(position_error_rightward)+1])
xlabel('Run number')
ylabel('Decoding error (cm)')
title({'Decoding error for all rightward trials', ['median: ' num2str(round(median(cell2mat(position_error_rightward))*100, 1)) ' cm']})

% Plot position errors per bin
for i = 1:length(rightward_t_start)
    len = length(true_binned_pos_rightward{i});
    [~, idx] = min(abs(repmat(true_binned_pos_rightward{i}, 66, 1) - repmat(linearised_bin_centres', 1, len)));
    true_binned_pos_binned_rightward{i} = linearised_bin_centres(idx);
end

% Find the decoding error in every linearised spatial bin
for i = 1:length(linearised_bin_centres)
    
    bin_error_rightward{i} = cell2mat(cellfun(@(x, y) x(y==linearised_bin_centres(i)), position_error_rightward, true_binned_pos_binned_rightward, 'UniformOutput', false));
    median_bin_error_rightward(i) = median(cell2mat(cellfun(@(x, y) x(y==linearised_bin_centres(i)), position_error_rightward, true_binned_pos_binned_rightward, 'UniformOutput', false)));
end
 
% Plot
figure; hold on

for i = 1:length(linearised_bin_centres)
    if ~isempty(bin_error{i})
        xjitter = linspace(linearised_bin_centres(i)-0.01, linearised_bin_centres(i)+0.01, length(bin_error{i}));
        scatter(xjitter, bin_error{i}*100, 'k.', 'MarkerEdgeAlpha', 0.1)
    end
end

plot(linearised_bin_centres, median_bin_error*100, 'b', 'LineWidth', 2)
xlabel('Linearised position (m)'); ylabel('Decoding error (cm)')
title('Distribution of decoding error across maze (all rightward runs)')

%% ===== Train decoder on all trials ===== %%

% Train decoder: get firing rate per bin
opts.smooth_kernel = 3;
opts.plot = false;
opts.adaptive_binning = true;
opts.spatial_bin_size = binsize;
[fr_per_bin_rightward, adaptive_bin_info_rightward] = trainDecoder(spiketimes, linearised_pos, linearised_bin_centres, Timestamps, rightward_t_start, rightward_t_end, opts);
[fr_per_bin_lefttward, adaptive_bin_info_leftward] = trainDecoder(spiketimes, linearised_pos, linearised_bin_centres, Timestamps, leftward_t_start, leftward_t_end, opts);

%% ==== REPLAY DECTECTION METHOD 1: decode activity during population bursts by fitting linear trajectories ==== %%

opts.binsize = 0.01;
opts.stepsize = 0.002;
opts.smooth_kernel = 3;
opts.plot = 'none'; % 'none', 'best_only', or 'as_you_go'
opts.nSh = 100;


gof_cell_shuffled_leftward = {};
gof_isi_shuffled_leftward = {};
gof_ratemap_shuffled_leftward = {};
gof_leftward = NaN(1, length(burst_start));
best_v_leftward = NaN(1, length(burst_start));
best_c_leftward = NaN(1, length(burst_start));

gof_cell_shuffled_rightward = {};
gof_isi_shuffled_rightward = {};
gof_ratemap_shuffled_rightward = {};
gof_rightward = NaN(1, length(burst_start));
best_v_rightward = NaN(1, length(burst_start));
best_c_rightward = NaN(1, length(burst_start));

opts.shuffle = [];

fprintf('%s: beginning decoding of rightward runs\n', datetime('now'))

parfor iBurst = 1:length(burst_start)
    
        [pxn, gof_leftward(iBurst), best_v_leftward(iBurst), best_c_leftward(iBurst), n] = decode(fr_per_bin_leftward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start(iBurst), burst_end(iBurst), opts);

        % Shuffle cell identities
        tmp = struct;
        tmp.binsize = 0.01; tmp.stepsize = 0.02; tmp.smooth_kernal = 3; tmp.plot = 'none'; tmp.nSh = 100; tmp.adaptive_binning = true; tmp.spatial_bin_size = binsize;
        tmp.shuffle = 'cell_identities';
        [gof_cell_shuffled_rightward{iBurst}] = decodeSignificance(fr_per_bin_rightward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start(iBurst), burst_end(iBurst), best_c_rightward(iBurst), best_v_rightward(iBurst), tmp);

        % Shuffle interspike intervals
        tmp.shuffle = 'interspike_intervals';
        [gof_isi_shuffled_rightward{iBurst}] = decodeSignificance(fr_per_bin_rightward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start(iBurst), burst_end(iBurst), best_c_rightward(iBurst), best_v_rightward(iBurst), tmp);

        % Circularly shift rate maps
        tmp.shuffle = 'circular_shufts';
        [gof_ratemap_shuffled_rightward{iBurst}] = decodeSignificance(fr_per_bin_rightward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start(iBurst), burst_end(iBurst), best_c_rightward(iBurst), best_v_rightward(iBurst), tmp);

        fprintf('%s: burst %i done\n', datetime('now'), iBurst)
    
end

fprintf('%s: beginning decoding of rightward runs\n', datetime('now'))

parfor iBurst = 1:length(burst_start)
    
        [pxn, gof_leftward(iBurst), best_v_leftward(iBurst), best_c_leftward(iBurst), n] = decode(fr_per_bin_leftward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start(iBurst), burst_end(iBurst), opts);

        % Shuffle cell identities
        tmp = struct;
        tmp.binsize = 0.01; tmp.stepsize = 0.02; tmp.smooth_kernal = 3; tmp.plot = 'none'; tmp.nSh = 100; tmp.adaptive_binning = true; tmp.spatial_bin_size = binsize;
        tmp.shuffle = 'cell_identities';
        [gof_cell_shuffled_leftward{iBurst}] = decodeSignificance(fr_per_bin_leftward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start(iBurst), burst_end(iBurst), best_c_rightward(iBurst), best_v_leftward(iBurst), tmp);

        % Shuffle interspike intervals
        tmp.shuffle = 'interspike_intervals';
        [gof_isi_shuffled_leftward{iBurst}] = decodeSignificance(fr_per_bin_leftward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start(iBurst), burst_end(iBurst), best_c_leftward(iBurst), best_v_leftward(iBurst), tmp);

        % Circularly shift rate maps
        tmp.shuffle = 'circular_shufts';
        [gof_ratemap_shuffled_leftward{iBurst}] = decodeSignificance(fr_per_bin_leftward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start(iBurst), burst_end(iBurst), best_c_leftward(iBurst), best_v_leftward(iBurst), tmp);

        fprintf('%s: burst %i done\n', datetime('now'), iBurst)
    
end

save([opts.save_name '_goodness_of_fit_rightward.mat'], 'gof_rightward', 'best_v_rightward', 'best_c_rightward', 'gof_cell_shuffled_rightward', 'gof_isi_shuffled_rightward', 'gof_ratemap_shuffled_rightward',...
    'gof_leftward', 'best_v_leftward', 'best_c_leftward', 'gof_cell_shuffled_leftward', 'gof_isi_shuffled_leftward', 'gof_ratemap_shuffled_leftward')

% Find significant replay
[sig_replay_rightward, loosely_sig_replay_rightward] = findSignificantReplays(gof_rightward, gof_cell_shuffled_rightward, gof_isi_shuffled_rightward, fr_per_bin_rightward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts);
disp([num2str(length(sig_replay_rightward)) ' significant replay events and ' num2str(length(loosely_sig_replay_rightward)) ' loosely significant replay events of rightward runs found'])

[sig_replay_leftward, loosely_sig_replay_leftward] = findSignificantReplays(gof_leftward, gof_cell_shuffled_leftward, gof_isi_shuffled_leftward, fr_per_bin_leftward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts);
disp([num2str(length(sig_replay_leftward)) ' significant replay events and ' num2str(length(loosely_sig_replay_leftward)) ' loosely significant replay events of leftward runs found'])

% Plot bursts with significant replay
plotSignificantReplay(sig_replay_rightward, fr_per_bin_rightward, best_v_rightward, best_c_rightward, gof_rightward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts)
plotSignificantReplay(sig_replay_leftward, fr_per_bin_leftward, best_v_leftward, best_c_leftward, gof_leftward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts)

% Plot when significant replays occur
plotReplayEvolution(Timestamps, rightward_t_start, rightward_t_end, Timestamps_q, burst_start, best_v, sig_replay_rightward, loosely_sig_replay_rightward, 'rightward');
plotReplayEvolution(Timestamps, leftward_t_start, leftward_t_end, Timestamps_q, burst_start, best_v, sig_replay_leftward, loosely_sig_replay_leftward, 'leftward');

%% ==== REPLAY DECTECTION METHOD 2: decode activity during population bursts using absolute weighted correlation ==== %%

% Decode population bursts
opts.binsize = 0.01;
opts.stepsize = 0.002;
opts.smooth_kernel = 3;
opts.plot = 'sig_only'; % 'none' or 'sig_only'
opts.nSh = 100;
opts.message = 'Decoding rightward runs';
[decoded_duration_rightward, decoded_arm_coverage_rightward, decoded_weighted_corr_rightward, decoded_replay_rightward] = decodeWeightedCorr(fr_per_bin_rightward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts);
opts.message = 'Decoding leftward runs';
[decoded_duration_leftward, decoded_arm_coverage_leftward, decoded_weighted_corr_leftward, decoded_replay_leftward] = decodeWeightedCorr(fr_per_bin_leftward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts);
fprintf('Found %i replays of rightward runs and %i replays of leftward runs\n', sum(decoded_replay_rightward), sum(decoded_replay_leftward))

% Shuffle
opts.shuffle = 'interspike_interval';
opts.nSh = 1000;
opts.message = 'Permuting';
[decoded_duration_sh, decoded_arm_coverage_sh, decoded_weighted_corr_sh, decoded_replay_sh, decoded_replay_combined_sh] = permutationWeightedCorr(fr_per_bin_rightward, fr_per_bin_leftward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts);

% Find p-values for cumulative replays
cum_replay_rightward = cumsum(decoded_replay_rightward);
cum_sh_rightward = cumsum(decoded_replay_sh(1:opts.nSh, :), 2);
pvalues_rightward = mean(cum_replay_rightward <= cum_sh_rightward);

cum_replay_leftward = cumsum(decoded_replay_leftward);
cum_sh_leftward = cumsum(decoded_replay_sh(1+opts.nSh:end, :), 2);
pvalues_leftward = mean(cum_replay_leftward <= cum_sh_leftward);

cum_replay_combined = cumsum(sum([decoded_replay_leftward; decoded_replay_rightward]));
cum_sh_combined = cumsum(decoded_replay_combined_sh, 2);
pvalues_combined = mean(cum_replay_combined <= cum_sh_combined);

% Plot cumulative replays
plotCumulativeReplays([decoded_replay_rightward + decoded_replay_leftward], decoded_replay_combined_sh, pvalues_combined, burst_start);

% Plot joint replays
joint = find(decoded_replay_rightward & decoded_replay_leftward);
plotJointReplays(fr_per_bin_leftward, fr_per_bin_rightward, spiketimes, linearised_bin_centres, Timestamps_q, burst_start(joint), burst_end(joint), opts)

% Save
save([opts.save_name '_weightedCorrelationReplay.mat'], 'decoded_duration_rightward', 'decoded_arm_coverage_rightward', 'decoded_weighted_corr_rightward', 'decoded_replay_rightward', 'decoded_duration_leftward', 'decoded_arm_coverage_leftward', 'decoded_weighted_corr_leftward', 'decoded_replay_leftward', 'decoded_duration_sh', 'decoded_arm_coverage_sh', 'decoded_weighted_corr_sh', 'decoded_replay_sh')
