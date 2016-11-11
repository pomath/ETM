function [st_series,etm] = adjust_polyhedra(st_series,poly,st_info,stab_sites)

    % remove stations with less than 50 epochs and those with < 1 yr of
    % data
    rm_stn = [];
    for i = 1:size(st_series,2)
        if size(st_series(i).epochs,2) < 50 | ...
                (max(st_series(2).epochs) - min(st_series(2).epochs)) < 1
            rm_stn = [rm_stn; st_series(i).stnm];
        end
    end

    for i = 1:size(rm_stn,1)
        [poly,st_series] = remove_stn(rm_stn(i,:), poly,st_series);
    end

    % find if there are any stab sites not present in the dataset
    stab_sites_i = zeros(size(stab_sites,1),1);
    for i = 1:length(stab_sites)
        index = structfind(st_series,'stnm',stab_sites{i});
        if ~isempty(index)
            stab_sites_i(i) = index;
        end
    end
    stab_sites_i = sortrows(stab_sites_i);
    stab_sites_i(stab_sites_i == 0) = [];

    % fit etms and filter large outliers
    % primera corrida de auto_fit_xyz busca datos que estén por afuera de 3
    % sigma y los descarta para no deformar mucho la red durante el ajuste
    [etm, index, weights] = auto_fit_xyz(st_series,st_info);

    % actualizo las series de tiempo sacando las colas estadisticas
    st_series = update_st_series(st_series,index,weights);

    [st_series, etm] = helmert_stacking(st_series, poly, st_info, etm, 4, stab_sites_i);

    for i = 1:4
        st_series_o = remove_cmm(st_series,poly,etm);
        [etm, index, weights] = auto_fit_xyz(st_series_o,st_info);
    end

    plot_etms(st_series, etm)

end