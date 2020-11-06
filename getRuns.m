
function [leftward_t_start, leftward_t_end, rightward_t_start, rightward_t_end] = getRuns(linearised_pos, speed, speed_thresh, Timestamps, Xpos, opts)

    % Threshold must be in same units as position (metres)

    leftward_t_start = [];
    leftward_t_end = [];
    rightward_t_start = [];
    rightward_t_end = [];
    speed(isnan(speed)) = 0;
    t = 1;
    
    % Remove any speed values where there is no matching linearised
    % position
    speed(isnan(linearised_pos(1:end-1))) = NaN;

    % Cycle through speed vector to find individual runs where speed is
    % above 50 cm/s
    while t < length(speed)

        % Find the next time that speed exceeds 0.3 m/s
        t_ = find(speed(t:end) > 0.3, 1) + t - 1;

        % Find start of run: last time speed was below threshold
        t_start = find(speed(t:t_) < speed_thresh | isnan(speed(t:t_)), 1, 'last') + t - 1;
        if isnan(speed(t_start)) t_start = t_start + 1; end

        % Find end of run: next time speed is below threshold (otherwise end of session)
        try t_end = find(speed(t_:end) < speed_thresh | isnan(speed(t_:end)), 1, 'first') + t_ - 1;
        catch t_end = length(speed); end
        if isnan(speed(t_end)) t_end = t_end - 1; end
        
        % Check that it covers at least 20% of the maze
        if range(Xpos(t_start:t_end)) / range(Xpos) < 0.2
            t = t_end + 1;
            continue
        end

        % Classify as leftward or rightward
        if linearised_pos(t_start) < linearised_pos(t_end)
            rightward_t_start = [rightward_t_start; t_start];
            rightward_t_end = [rightward_t_end; t_end];

        elseif linearised_pos (t_start) > linearised_pos(t_end)
            leftward_t_start = [leftward_t_start; t_start];
            leftward_t_end = [leftward_t_end; t_end];
        end

        t = t_end + 1;

    end
    
    % Plot
    if isfield(opts, 'plot') && opts.plot == true
        
        figure; hold on
        
        % Plot X-position
        plot(Timestamps, Xpos, '.k')
        axis tight;
        
        % Plot times of runs
        for iRun = 1:length(rightward_t_start)
            p = patch([Timestamps(rightward_t_start(iRun)) Timestamps(rightward_t_start(iRun)) Timestamps(rightward_t_end(iRun)) Timestamps(rightward_t_end(iRun))],...
                [min(Xpos) max(Xpos) max(Xpos) min(Xpos)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
            if iRun > 1 set(p, 'HandleVisibility', 'off'); end
        end
        for iRun = 1:length(leftward_t_start)
            p = patch([Timestamps(leftward_t_start(iRun)) Timestamps(leftward_t_start(iRun)) Timestamps(leftward_t_end(iRun)) Timestamps(leftward_t_end(iRun))],...
                [min(Xpos) max(Xpos) max(Xpos) min(Xpos)], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
            if iRun > 1 set(p, 'HandleVisibility', 'off'); end
        end
        
        title('Runs detected')
        l = legend('Animal''s X-position', 'Detected rightward runs', 'Detected leftward runs', 'Location', 'southoutside'); set(l, 'Box', 'off');
    
    end
end