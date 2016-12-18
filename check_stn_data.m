function r = check_stn_data(t,level)
    % stations with less than 50 epochs and those with < 1 yr of data
    % and if the actual data is < 10% of the estimated data size, 
    % should not have an ETM adjustment
    % this avoids wild excursions of periodic terms due to impossibility of
    % adjusting the correct truncated fourier series
    ts = min(t);
    te = max(t);
    % estimate the number of epochs for a full time series
    epochs = (te-ts)*365;
    % find the time gaps between epochs
    dt = diff(t);
    if length(t) < 50 | (te - ts) < 1 | level*epochs > sum(dt<0.5)
        r = false;
    else
        r = true;
    end
end