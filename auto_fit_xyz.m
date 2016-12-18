function [etm, oindex, weights] = auto_fit_xyz(st_series,st_info,start_at)

    % initialize the cell arrays to save the etm data
    etm = cell(size(st_series,2),2);
    % initialize the cell array for the outlier index information
    oindex = cell(size(st_series,2),1);
    % weights for the GPS data
    weights = cell(size(st_series,2),1);
    
    % station counter
    if nargin == 3
        if size(start_at,1) == 1
            j = (find(strcmp({st_series.stnm}, start_at)==1):length(st_series))';
        else
            j = sort(start_at);
        end
    else
        j = (1:length(st_series))';
    end

    for i = j'

        stnm = st_series(i).stnm;

        % read in antenna changes from station info and co-seismic jumps
        check_files(stnm, st_info, st_series(i).lat, st_series(i).lon);
        
        % get the epochs for the current station
        t = (st_series(i).epochs)';
        
        % run a check to see if you can fit something to the dataset
        % if false, this stations does not have the min requirements to fit
        % an ETM, skip it
        
        if check_stn_data(t,0.2)
        
            % get the latitude of current station
            % this is only necessary for the TIDES code, no need for the regular GPS ETM
            lat = st_series(i).lat;
            lon = st_series(i).lon;

            % save the vector with the observations (for simplicity)
            % all observations necessary for the fit
            Lx = (st_series(i).x)';
            Ly = (st_series(i).y)';
            Lz = (st_series(i).z)';

            % after taking to Mike, we don't need to use the formal errors from GAMIT, since they don't mean anything
            % so use unit weights and then update during iteration using sigma zero information
            px = ones(size(st_series(i).px));
            py = ones(size(st_series(i).px));
            pz = ones(size(st_series(i).px));

            % load the heavy-side functions and periodic terms
            [A,Ha,constrains_h,constrains_s,~] = load_hsf3(stnm, t,true);

            if cond(A) > 1e3
                disp(['  >> Unstable system of equations for station: ' stnm '. It has been stabilized, but make sure to check the data and results.'])
            end

            for comp = 1:3
                switch comp
                    case 1
                        p = px; L = Lx;
                    case 2
                        p = py; L = Ly;
                    case 3
                        p = pz; L = Lz;
                end

                % weights
                P = diag(1./p.^2);

                % print current component
                fprintf('[%i]',comp);

                [C, S, So, V, ~, ~, cst_pass, p, index] = adjust_lsq(A,P,L,constrains_h,constrains_s);

                if cst_pass
                    cst_pass = 'PASS';
                else
                    cst_pass = 'NOT PASS';
                end

                % update the information
                switch comp
                    case 1
                        Cx = C; Sx = S; sox = So; Vx = V; Px = p; xresult = cst_pass; index_x = index;
                    case 2
                        Cy = C; Sy = S; soy = So; Vy = V; Py = p; yresult = cst_pass; index_y = index;
                    case 3
                        Cz = C; Sz = S; soz = So; Vz = V; Pz = p; zresult = cst_pass; index_z = index;
                end
            end

            fprintf('\n');
            disp([stnm ' sigma0_north = ' sprintf('%f', sox) ' sigma0_east = ' sprintf('%f', soy) ' sigma0_up = ' sprintf('%f', soz) ' ' xresult ' ' yresult ' ' zresult])
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % pack everything in the cell array
            oindex(i) = num2cell(index_x & index_y & index_z ,1);

            weights(i) = num2cell([Px'; Py'; Pz'],[1 2]);
            etm(i,1) = num2cell([Cx'; Cy'; Cz'],[1 2]);
            etm(i,2) = num2cell(Ha,[1 2]);
        end
        % save_etm_params(adjcmp,min(t),max(t),[Cx Cy Cz],Ha,stnm)
    end
end
