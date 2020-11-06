
function [Xpos, Ypos, Timestamps] = processPositionData(Xpos, Ypos, Timestamps, MazeEpoch, start_index, end_index)

    % Restrict to start_index and end_index given by user
    Xpos = Xpos(start_index:end_index);
    Ypos = Ypos(start_index:end_index);
    Timestamps = Timestamps(start_index:end_index);
    
    % Restrict to the period defined as the maze epoch
    Xpos = Xpos(Timestamps > MazeEpoch(1) & Timestamps < MazeEpoch(2));
    Ypos = Ypos(Timestamps > MazeEpoch(1) & Timestamps < MazeEpoch(2));
    Timestamps = Timestamps(Timestamps > MazeEpoch(1) & Timestamps < MazeEpoch(2));
    
    % Interpolate position data to get a position for every milisecond
    Timestamps = round(Timestamps, 3);
    Xpos = interp1(Timestamps,Xpos,Timestamps(1):0.001:Timestamps(end));
    Ypos = interp1(Timestamps,Ypos,Timestamps(1):0.001:Timestamps(end));
    Timestamps = Timestamps(1):0.001:Timestamps(end);

end