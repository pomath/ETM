function poly = create_poly_struct(st_series)
% create the polyhedra structure using the st_series struct

    % build the polyhedra structure
    poly.epochs = sortrows(unique([st_series.epochs]))';
    % make a vector for coordinates (rows time, cols station
    poly.x = nan(size(poly.epochs,1),size(st_series,2));
    poly.y = nan(size(poly.epochs,1),size(st_series,2));
    poly.z = nan(size(poly.epochs,1),size(st_series,2));
    poly.px = nan(size(poly.epochs,1),size(st_series,2));
    poly.py = nan(size(poly.epochs,1),size(st_series,2));
    poly.pz = nan(size(poly.epochs,1),size(st_series,2));
    % use the poly.epochs vector to sort the time series into polyhedra
    for i = 1:size(st_series,2)
        c = ismember(poly.epochs,st_series(i).epochs);
        poly.x(c,i) = st_series(i).x;
        poly.y(c,i) = st_series(i).y;
        poly.z(c,i) = st_series(i).z;
        poly.px(c,i) = st_series(i).px;
        poly.py(c,i) = st_series(i).py;
        poly.pz(c,i) = st_series(i).pz;
    end
end

