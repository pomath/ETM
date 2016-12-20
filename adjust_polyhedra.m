function [st_series,etm] = adjust_polyhedra(st_series,poly,st_info,stab_sites_i)

    % fit etms and filter large outliers
    % do it only for the stab sites
    [etm, index, weights] = auto_fit_xyz(st_series,st_info,stab_sites_i);

    % actualizo las series de tiempo sacando las colas estadisticas
    st_series = update_st_series(st_series,index,weights,stab_sites_i);

    [st_series, etm, ~] = helmert_stacking(st_series, poly, st_info, etm, 4, stab_sites_i);

    % remove the common mode
    st_series = remove_cmm(st_series,etm);

    % readjust the ETM without the common mode
    [etm, index, weights] = auto_fit_xyz(st_series,st_info);
    st_series = update_st_series(st_series,index,weights);
    
    % plot the stab sites only
    plot_etms(st_series, etm, 'etm_stab', stab_sites_i)
    
    save('../ETM_files/st_series/st_series.mat','st_series');
    save('../ETM_files/st_series/etm.mat','etm');
    save('../ETM_files/st_series/stab_sites_i.mat','stab_sites_i');
    save('../ETM_files/st_series/st_info.mat','st_info');
end