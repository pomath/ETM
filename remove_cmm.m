function st_series_o = remove_cmm(st_series_r,poly,etm)
    % remove common mode
    % each station has an XYZ coordinate and S1,C1,S2,C2 parameter
    % want to find the rotation and translation that minimizes S1,C1,S2,C2

    st_series_o = st_series_r;
    
    X = []; Y = []; Z = [];
    % loop through the stations and get one coordinate for each one
    for i = 1:size(st_series_r,2)
        X(1,i) = mean(st_series_r(i).x);
        Y(1,i) = mean(st_series_r(i).y);
        Z(1,i) = mean(st_series_r(i).z);
    end

    % build the design matrix
    Ax = [zeros(size(poly.x,2),1) Z' -Y'];
    Ay = [-Z' zeros(size(poly.x,2),1) X'];
    Az = [Y' -X' zeros(size(poly.x,2),1)];

    A = [Ax; Ay; Az];

    for i = 1:4
        % L has the parameters S1,C1,S2,C2
        L = [osc_param(i,:)'; osc_param(i+4,:)'; osc_param(i+8,:)'];
        x = (A'*A)\(A'*L);

        % remove from each time series
        % the translation part has to be directly removed from the time series
        switch i
            case 1
                t = sin(2*pi.*poly.epochs);
            case 2
                t = cos(2*pi.*poly.epochs);
            case 3
                t = sin(4*pi.*poly.epochs);
            case 4
                t = cos(4*pi.*poly.epochs);
        end

        for j = 1:size(st_series_r,2)
            st_series_o(j).x = st_series_o(j).x - (t(ismember(poly.epochs,st_series_r(j).epochs)).*mean(osc_param(i,:)))';
            st_series_o(j).y = st_series_o(j).y - (t(ismember(poly.epochs,st_series_r(j).epochs)).*mean(osc_param(i+4,:)))';
            st_series_o(j).z = st_series_o(j).z - (t(ismember(poly.epochs,st_series_r(j).epochs)).*mean(osc_param(i+8,:)))';

            % get the amplitude of the rotational component
            a = x(1).*X(j);
            b = x(2).*Y(j);
            c = x(3).*Z(j);

            st_series_o(j).x = st_series_o(j).x - (t(ismember(poly.epochs,st_series_r(j).epochs)).*a)';
            st_series_o(j).y = st_series_o(j).y - (t(ismember(poly.epochs,st_series_r(j).epochs)).*b)';
            st_series_o(j).z = st_series_o(j).z - (t(ismember(poly.epochs,st_series_r(j).epochs)).*c)';
        end
    end
end

