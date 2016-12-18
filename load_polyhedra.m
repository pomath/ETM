function [poly,st_series] = load_polyhedra(folder)
    % function to load the GPS polyhedra and GPS time series
    %   it returns 2 structures:
    %       poly: the structure with the epochs and XYZ coordinates
    %           where the row index is the date and the cols are stations 
    %       st_series: classic rtvel structure with the time series
    %
    % initialize the struct for the time series
    % this data comes from make_data script
    % DDG: 12/06/2016
    % esta versión es diferente de la clásica que utilizo para los ETM del
    % IGN porque tiene un agregado para leer las componentes de marea de
    % cada estación
       
    st_series(1,1).stnm = [];
    st_series(1,1).epochs = [];
    st_series(1,1).n = [];
    st_series(1,1).e = [];
    st_series(1,1).d = [];
    st_series(1,1).x = [];
    st_series(1,1).y = [];
    st_series(1,1).z = [];
    st_series(1,1).pn = [];
    st_series(1,1).pe = [];
    st_series(1,1).pd = [];
    st_series(1,1).px = [];
    st_series(1,1).py = [];
    st_series(1,1).pz = [];
    st_series(1,1).lat = [];
    st_series(1,1).lon = [];
    
    eesum = [];
    
    % load the time series files
    files=dir([folder '/*.txt']);

    ss=size(st_series,2);
    for i = 1:size(files,1)
        % los archivos de las series de tiempo son construidos por
        % make_data(v) v = version
        fileID = fopen([folder '/' files(i).name],'r');

        disp(files(i).name)
        
        % fgetl(fileID);
        stations = textscan(fileID, '%s %f %f %f %f %f %f %f %f %f %f %f %f %f');
        
        if size(stations{9}',2) == 0
            disp 'stop'
        end
       
        st_series(1,i).stnm = files(i).name(1:4);
        st_series(1,i).epochs = stations{2}';
        % all indexes used (no knowledge of outliers)
        st_series(1,i).index = true(size(stations{9},1),1);
        
        st_series(1,i).n = stations{3}';
        st_series(1,i).e = stations{4}';
        st_series(1,i).d = stations{5}';
        st_series(1,i).pn = stations{6}';
        st_series(1,i).pe = stations{7}';
        st_series(1,i).pd = stations{8}';
        
        st_series(1,i).x = stations{9}';
        st_series(1,i).y = stations{10}';
        st_series(1,i).z = stations{11}';
        st_series(1,i).px = stations{12}';
        st_series(1,i).py = stations{13}';
        st_series(1,i).pz = stations{14}';

        x = stations{9}';
        y = stations{10}';
        z = stations{11}';

        [lat, lon, ~] = ecef2lla(x(1), y(1), z(1));
        st_series(1,i).lat = lat*180/pi;
        if lon*180/pi > 180
            lon=lon*180/pi-360;
        else
            lon=lon*180/pi;
        end
        st_series(1,i).lon = lon;
        
        if ~isempty(x == 0) | ~isempty(y == 0) | ~isempty(z == 0)
            % remove nonsense values
            st_series(1,i).x(x == 0 | y == 0 | z == 0) = [];
            st_series(1,i).y(x == 0 | y == 0 | z == 0) = [];
            st_series(1,i).z(x == 0 | y == 0 | z == 0) = [];
            st_series(1,i).px(x == 0 | y == 0 | z == 0) = [];
            st_series(1,i).py(x == 0 | y == 0 | z == 0) = [];
            st_series(1,i).pz(x == 0 | y == 0 | z == 0) = [];
            st_series(1,i).pn(x == 0 | y == 0 | z == 0) = [];
            st_series(1,i).pe(x == 0 | y == 0 | z == 0) = [];
            st_series(1,i).pd(x == 0 | y == 0 | z == 0) = [];
            st_series(1,i).n(x == 0 | y == 0 | z == 0) = [];
            st_series(1,i).e(x == 0 | y == 0 | z == 0) = [];
            st_series(1,i).d(x == 0 | y == 0 | z == 0) = [];
            st_series(1,i).epochs(x == 0 | y == 0 | z == 0) = [];
            st_series(1,i).index(x == 0 | y == 0 | z == 0) = [];
        end
    end

    % la estructura de los poliedros es un poco diferente, dado que aunque
    % no haya datos para una epoca, los vectores deben tener todos el mismo
    % tamaño.
    poly = create_poly_struct(st_series);
    
end