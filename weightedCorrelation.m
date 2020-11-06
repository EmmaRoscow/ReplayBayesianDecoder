
function [weighted_corr] = weightedCorrelation(pxn, fr_per_bin, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts)

    % Get posterior probabilities if not supplied
    if isempty(pxn)
        pxn = decode(fr_per_bin, spiketimes, linearised_bin_centres, Timestamps_q, burst_start, burst_end, opts);
    end
    
    % Find weighted correlation
    x = repmat([1:size(pxn, 1)]', 1, size(pxn, 2));       % Array of time bins
    y = repmat([1:size(pxn, 2)], size(pxn, 1), 1);        % Array of position bins
    m_xw = sum(sum(pxn .* x)) / sum(sum(pxn));
    m_yw = sum(sum(pxn .* y)) / sum(sum(pxn));
    
    cov_xyw = sum(sum(pxn .* ( x - m_xw ) .* ( y - m_yw))) / sum(sum(pxn));
    cov_xxw = sum(sum(pxn .* ( x - m_xw ) .* ( x - m_xw))) / sum(sum(pxn));
    cov_yyw = sum(sum(pxn .* ( y - m_yw ) .* ( y - m_yw))) / sum(sum(pxn));
    
    weighted_corr = cov_xyw / sqrt( cov_xxw * cov_yyw );
    
end