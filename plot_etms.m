function plot_etms(st_series, etm, start_at)

    % station counter
    if nargin == 3
        i = structfind(st_series,'stnm',start_at);
    else
        i = 1;
    end

    while i <= size(st_series,2)

        % station name
        stnm = st_series(i).stnm;
        index = st_series(i).index;
        
        % get the epochs for the current station
        t = (st_series(i).epochs)';
        % get the latitude of current station
        lat = st_series(i).lat;
        lon = st_series(i).lon;
        
        Lx = (st_series(i).x);
        Ly = (st_series(i).y);
        Lz = (st_series(i).z);
        
        % ETM parameters
        Ha = etm{i,2};
        C = etm{i,1};
        Cx = C(1,:);
        Cy = C(2,:);
        Cz = C(3,:);
        
        [Cn,Ce,Cu] = ct2lg(Cx',Cy',Cz',lat*pi/180,lon*pi/180);
        [Ln,Le,Lu] = ct2lg(Lx',Ly',Lz',lat*pi/180,lon*pi/180);
        
        plot_station(stnm, i, index, lat, lon, t, Ha, Cn, Ce, Cu, Ln, Le, Lu)
        
        i = i + 1;
    end

end

