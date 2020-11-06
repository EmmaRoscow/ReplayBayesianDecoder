
function [linearised_pos, linearised_bin_centres, pos_bins] = linearisePositions(Xpos, Ypos, binsize, mazetype)

    switch mazetype

        case '1.6m Linear Maze'

            Xbin = floor(min(Xpos)/binsize)*binsize : binsize : ceil(max(Xpos)/binsize)*binsize;

            % Snap every position to a bin centre
            Xpos_binned = double.empty(length(Xpos), 0);
            for i = 1:length(Xpos)
                % Keep NaNs as NaNs
                if isnan(Xpos(i))
                    Xpos_binned(i) = NaN;
                else
                    % Snap to Xbin
                    distance = abs(Xpos(i) - Xbin);
                    [~, closest_bin(i)] = min(distance);
                    Xpos_binned(i) = Xbin(closest_bin(i));
                end

            end
            
            linearised_pos = Xpos_binned;
            linearised_bin_centres = Xbin;
            pos_bins = closest_bin;
            
        otherwise
            error(['Don''t know how to process this maze type: ' mazetype])
    
    end
    
end