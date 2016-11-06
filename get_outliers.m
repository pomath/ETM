function out = get_outliers(st_series,epoch)
    out = [];
    for i = 1:size(st_series,2)
        if ismember(epoch,st_series(i).epochs)
            out(i,1) = st_series(i).index(st_series(i).epochs == epoch);
        else
            out(i,1) = 0;
        end
    end
end

