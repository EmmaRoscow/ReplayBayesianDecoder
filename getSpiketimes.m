
function spiketimes = getSpiketimes(SpikeTimes, SpikeIDs, PyrIDs)

    % SpikeTimes is a vector of spiketimes in seconds (for all units)

    % SpikeIDs is a vector with same length as SpikeTimes, indicating the ID of
    % the cell which emitted the corresponding spike

    % PyrIDs is a vector of unique IDs for pyramidal cells, generally much
    % shorter than the other two vectors

    % Initialise cell array for spiketimes (one cell per pyramidal cell)
    spiketimes = cell(length(PyrIDs), 1);

    % For each pyramidal cell, find all spikes which correspond to it
    for i = 1:length(PyrIDs)

        iPyr = PyrIDs(i);
        spiketimes{i} = SpikeTimes(SpikeIDs == iPyr);

    end

end