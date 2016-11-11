function out = get_outliers_stab(st_series,epoch,stab_sites)
    out = [];
    for i = 1:size(st_series,2)
        if ismember(epoch,st_series(i).epochs) & ismember(i,stab_sites)
            out(i,1) = st_series(i).index(st_series(i).epochs == epoch);
        else
            out(i,1) = 0;
        end
    end
end

