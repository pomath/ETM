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
        
        % do not load the Ha from file; use the Ha from the real data
        % adjustment
        % do not send ts in full. Will cause problems with heavy side
        % functions from log decay that occured before there was data in
        % the time series only calculate A with data within limits and fill
        % the rest with NaNs
        t = st_series(i).epochs';
    
        [A,~] = load_hsf3(stnm, t, false,Ha);

        C = etm{i,1};
        Cx = C(1,:)';
        Cy = C(2,:)';
        Cz = C(3,:)';
        
        if size(A,2) ~=size(Cx,1)
            disp(['Warning: ' stnm ' data and A matrix size do not agree'])
        end
        
        tindex = ismember(ts,t);
        
        Xe(i,tindex) = A*Cx;
        Ye(i,tindex) = A*Cy;
        Ze(i,tindex) = A*Cz;
    end
end