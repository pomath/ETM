function [Xe,Ye,Ze] = get_coords_etm(st_series,etm,ts)
% esta funcion calcula las coordenadas de las estaciones utilizando los
% parámetros estimados del ETM
%    
    Xe = nan(size(st_series,2),size(ts,1));
    Ye = nan(size(st_series,2),size(ts,1));
    Ze = nan(size(st_series,2),size(ts,1));

    for i = 1:size(st_series,2)
        stnm = st_series(i).stnm;
        Ha = etm{i,2};
        lat = st_series(i).lat;
        
        % do not load the Ha from file; use the Ha from the real data
        % adjustment
        [A,~] = load_hsf3(stnm, ts, lat,Ha);

        C = etm{i,1};
        Cx = C(1,:)';
        Cy = C(2,:)';
        Cz = C(3,:)';
        
        Xe(i,:) = A*Cx;
        Ye(i,:) = A*Cy;
        Ze(i,:) = A*Cz;
    end
end