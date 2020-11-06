
function speed = getSpeed(Xpos, Ypos, Timestamps, kernel)

    % Smooth position data (depending on Matlab version)
    try
        Xpos_smoothed = smoothts(Xpos, 'b', kernel);
        Ypos_smoothed = smoothts(Ypos, 'b', kernel);
    catch
        Xpos_smoothed = smoothdata(Xpos, 'movmean', kernel);
        Ypos_smoothed = smoothdata(Ypos, 'movmean', kernel);
    end

    % Calculate speed at every timestamp: Euclidean distance between
    % consecutive timestamps divided by time difference
    for t = 1:length(Timestamps)-1
        speed(t) = sqrt( (Xpos_smoothed(t+1)-Xpos_smoothed(t)).^2 + (Ypos_smoothed(t+1)-Ypos_smoothed(t)).^2 ) / (Timestamps(t+1)-Timestamps(t));
    end

end