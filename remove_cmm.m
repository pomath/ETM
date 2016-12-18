function st_series_o = remove_cmm(st_series_r,etm)
    % remove common mode
    % each station has an XYZ coordinate and S1,C1,S2,C2 parameter
    % want to find the rotation and translation that minimizes S1,C1,S2,C2
    
    % rebuild the poly structure to eliminate possible epochs removed by
    % helmert_stacking
    poly = create_poly_struct(st_series_r);
    
    X = []; Y = []; Z = [];
    osc_x = []; osc_y = []; osc_z = [];
    % loop through the stations and get one coordinate for each one
    for i = 1:size(st_series_r,2)
        % use 50% of data presence in check_stn_data
        if ~isempty(etm{i,1}) & check_stn_data(st_series_r(i).epochs,0.5)
            osc_x = [osc_x; etm{i,1}(1,end-3:end)];
            osc_y = [osc_y; etm{i,1}(2,end-3:end)];
            osc_z = [osc_z; etm{i,1}(3,end-3:end)];
            index(i) = 1;
        else
            index(i) = 0;
        end
        % get a mean value for all stations, even if it doesn't have an etm
        X(1,i) = mean(st_series_r(i).x);
        Y(1,i) = mean(st_series_r(i).y);
        Z(1,i) = mean(st_series_r(i).z);
    end
    
    index = logical(index);
    
    n = size(X,2);
    
    % build the design matrix
    Ax = [zeros(n,1) -Z' Y' repmat([1 0 0],n,1)];
    Ay = [Z' zeros(n,1) -X' repmat([0 1 0],n,1)];
    Az = [-Y' X' zeros(n,1) repmat([0 0 1],n,1)];

    A = [Ax; Ay; Az];
    Ai = [Ax; Ay; Az];
    
    % remove the stations with no ETM
    A([index'; index'; index'] == 0,:) = [];
    
    rt = [poly.x poly.y poly.z];
    
    for i = 1:4
        
        % L has the parameters S1,S2,C1,C2
        L = [osc_x(:,i); osc_y(:,i); osc_z(:,i)];
        x = (A'*A)\(A'*L);

        % remove from each time series
        % the translation part has to be directly removed from the time series
        switch i
            case 1
                t = repmat(sin(2*pi.*poly.epochs),1,size(poly.x,2));
            case 2
                t = repmat(sin(4*pi.*poly.epochs),1,size(poly.x,2));
            case 3
                t = repmat(cos(2*pi.*poly.epochs),1,size(poly.x,2));
            case 4
                t = repmat(cos(4*pi.*poly.epochs),1,size(poly.x,2));
        end
        % remove the common mode for this component
        rt = rt - [t t t].*repmat((Ai*x)',size(t,1),1);
    end
       
    st_series_o = load_poly2series(st_series_r,rt,[]);
end

