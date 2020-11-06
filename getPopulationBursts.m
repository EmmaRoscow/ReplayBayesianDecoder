
function [burst_start, burst_end, quiescence_mua, Timestamps_q] = getPopulationBursts(Timestamps, spiketimes, speed)

    % Find when speed goes above or below 0.05 m/s
    speed_change = diff(speed < 0.05);
    
    % Find start of quiescence (when speed goes below 0.05 m/s) and end of
    % quiescence (when speed goes above)
    quiescence_start = find(speed_change == 1);
    quiescence_end = find(speed_change == -1);
    if quiescence_start(1) > quiescence_end(1)
        quiescence_end = quiescence_end(2:end);
    end
    if length(quiescence_end) < length(quiescence_start)
        quiescence_end = [quiescence_end length(speed_change)];
    end
    
    % Restrict to quiescence bouts which are long enough (at least 1 second)
    long_enough = find(Timestamps(quiescence_end) - Timestamps(quiescence_start) >= 1);
    quiescence_start = quiescence_start(long_enough);
    quiescence_end = quiescence_end(long_enough);
    
    % Find multi-unit activity in 1ms bins during quiescence bouts
    quiescence_mua = {};
    Timestamps_mua = {};

    f = waitbar(0, ['Calculating multi-unit activity...']); tic
    
    for i = 1:length(quiescence_start)
    
        % Get spikes in 10ms bins
        Timestamps_mua{i} = round(Timestamps(quiescence_start(i)), 3) : 0.01 : round(Timestamps(quiescence_end(i)), 3);
        bout_mua = [];
        for tt = round(Timestamps(quiescence_start(i)), 3) : 0.01 : round(Timestamps(quiescence_end(i)), 3)
            bout_mua = [bout_mua sum(spiketimes >= tt & spiketimes < tt+0.01)];
            waitbar((i-1)/length(quiescence_start), f, {'Calculating multi-unit activity...', [num2str(round(toc)) ' seconds elapsed']});
        end

        quiescence_mua{i} = bout_mua;

    end
    
    % Find when multi-unit activity exceeds threshold of mean + 3 standard deviations
    mean_mua = mean(cell2mat(quiescence_mua));
    std_mua = std(cell2mat(quiescence_mua));
        mua_change = diff(cell2mat(quiescence_mua) > mean_mua + 3*std_mua);
    burst_peak = find(mua_change == 1);

    burst_start = [];
    burst_end = [];

    all_quiescence_mua = cell2mat(quiescence_mua);

    % Waitbar to plot progress
    waitbar(0, f, 'Finding MUA population bursts...');

    % Find start and end of individual population bursts
    for i = 1:length(burst_peak)
        next_start = find(all_quiescence_mua(1:burst_peak(i)) < mean_mua, 1, 'last');
        if i == 1 || next_start > burst_end(end)
            burst_start = [burst_start next_start];
            burst_end = [burst_end find(all_quiescence_mua(burst_peak(i):end) < mean_mua, 1, 'first') + burst_peak(i) - 1];
        end
        waitbar(i / length(burst_peak), f, {'Finding MUA population bursts...', [num2str(round(toc)) ' seconds elapsed']});
    end
    
    % Record timestamps during quiescence
    Timestamps_q = cell2mat(Timestamps_mua);
    
    % Restrict to bursts which last 50 - 400 milliseconds
    burst_duration = Timestamps_q(burst_end) - Timestamps_q(burst_start);
    burst_start = burst_start(burst_duration >= 0.05 & burst_duration <= 0.4);
    burst_end = burst_end(burst_duration >= 0.05 & burst_duration <= 0.4);
    
end
