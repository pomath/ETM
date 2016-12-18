function [stab_sites,stab_sites_i] = load_stab_sites(st_series,file)
    fileID = fopen(file,'r');

    stab_sites = textscan(fileID,'%s');
    stab_sites = stab_sites{1};
    
    fclose(fileID);
    
    if size(stab_sites,1) == 1
        if strcmp(stab_sites,'ALL')
            stab_sites = {st_series.stnm}';
        end
    end
    
    % find if there are any stab sites not present in the dataset
    stab_sites_i = zeros(size(stab_sites,1),1);
    for i = 1:length(stab_sites)
        index = find(strcmp({st_series.stnm}, stab_sites{i})==1);
        if ~isempty(index)
            % check if this can be a stabilization site (need to be able to
            % have an ETM and > 50% of the data for its timespan)
            if check_stn_data(st_series(index).epochs,0.5)
                stab_sites_i(i) = index;
            end
        end
    end
    stab_sites_i = sortrows(stab_sites_i);
    stab_sites_i(stab_sites_i == 0) = [];
end