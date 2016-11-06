function st_series = load_daily()
    load '~/Dropbox/Geofisica/Dissertation/maule_model/gps_stations/maule_SA.mat'; ts_rtvel4 = st_series;

    st_series = [];
    st_lla = [];
    % select the data we are interested in
    for i = 1:size(ts_rtvel4,2)
        if (ts_rtvel4(1,i).lonlat(2) < -26 && ts_rtvel4(1,i).lonlat(2) > -45) && ...
                (ts_rtvel4(1,i).lonlat(1) < -55 && ts_rtvel4(1,i).lonlat(1) > -80) && ...
                size(ts_rtvel4(1,i).epochs,2) >= 360 && ...
                (min(ts_rtvel4(1,i).epochs) < 2010 && max(ts_rtvel4(1,i).epochs) > 2010.5)

            % this is data that we want, save it
            st_lla = [st_lla; ts_rtvel4(1,i).lonlat(2) ts_rtvel4(1,i).lonlat(1)];
            st_series = [st_series ts_rtvel4(1,i)];
        end
    end
end