function out = get_outliers_stab(st_series,poly,stab_sites)
% create a matrix with station-rows and column-epochs with the outlier
% information. If a station is not in the stab_sites vector, the column is
% all zeros

    out = false(size(st_series,2),size(poly.epochs,1));
    
    for i = 1:size(st_series,2)
        if ismember(i,stab_sites)
            c = ismember(poly.epochs,st_series(i).epochs);
            out(i,c) = st_series(i).index;
        end
    end
end

