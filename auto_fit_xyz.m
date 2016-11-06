function [etm, oindex, weights] = auto_fit_xyz(st_series,st_info,start_at)

    etm = cell(size(st_series,2),2);
    oindex = cell(size(st_series,2),1);
    weights = cell(size(st_series,2),1);
    
    % station counter
    if nargin == 3
        i = structfind(st_series,'stnm',start_at);
    else
        i = 1;
    end

    while i <= size(st_series,2)

        stnm = st_series(i).stnm;

        % configura los archivos con saltos (lee sismos y station_info)
        check_files(stnm, st_info, st_series(i).lat, st_series(i).lon);

        % be aware of any outliers
        index = st_series(i).index;
        
        % get the epochs for the current station
        t = (st_series(i).epochs)';
        
        % get the latitude of current station
        lat = st_series(i).lat;
        lon = st_series(i).lon;
        
        % save the vector with the observations (for simplicity)
        % all observations necessary for the fit
        Lx = (st_series(i).x)';
        Ly = (st_series(i).y)';
        Lz = (st_series(i).z)';

        px = st_series(i).px;
        py = st_series(i).py;
        pz = st_series(i).pz;
        
        % load the heavy-side functions and periodic terms
        [A,Ha] = load_hsf3(stnm, t, lat);
        
        for comp = 1:3
            switch comp
                case 1
                    p = px; L = Lx;
                case 2
                    p = py; L = Ly;
                case 3
                    p = pz; L = Lz;
            end

            % factor for downweight or upweight of P
            factor = 1;
            % factor to downweight outliers
            factor_outliers = 1;
            
            % weights
            P = diag(ones(1,size(p,2)));
            
            % goodness of fit test
            cst_pass = false;
            % print current component
            fprintf('[%i]',comp);
            
            for y=1:2
                % funcion de ajuste por MMCC personalizada. Además de
                % ajustar hace un test de bondad de ajuste y muestra si
                % pasa o no pasa. Se ejecuta iterativamante corrigiendo los
                % pesos hasta que pase el test.
                [C, S, So, V, r, dof, cst_pass] = adjust_lsq(A,P,L,index);
                
                % verify if the test de bondad was passed
                   
                % determine the outliers
                [~,outliers]=deleteoutliers(V, 0.01);
                % downweight them (x2)
                factor_outliers=ones(size(p));
                factor_outliers(outliers) = 10;
                % make the new P matrix
                P = diag(1./((ones(1,size(p,2)).*factor_outliers).^2));
            end

            switch comp
                case 1
                    Cx = C; Sx = S; sox = So; Vx = V; Px = p.*factor.*factor_outliers;
                case 2
                    Cy = C; Sy = S; soy = So; Vy = V; Py = p.*factor.*factor_outliers;
                case 3
                    Cz = C; Sz = S; soz = So; Vz = V; Pz = p.*factor.*factor_outliers;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % test de bondad de ajuste
        % para mostrar en matlab solamente
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        xx = sox.^2.*dof;
        xy = soy.^2.*dof;
        xz = soz.^2.*dof;
        % careful! This function returns the opposite value of alpha as on
        % Leick, page 143
        X1 = chi2inv(1-0.05/2,dof);
        X2 = chi2inv(0.05/2,dof);

        if xx < X2 || xx > X1; xresult = 'FAIL'; else xresult = 'PASS'; end
        if xy < X2 || xy > X1; yresult = 'FAIL'; else yresult = 'PASS'; end
        if xz < X2 || xz > X1; zresult = 'FAIL'; else zresult = 'PASS'; end

        fprintf('\n');
        disp([stnm ' sigma0_north = ' num2str(sox) ' sigma0_east = ' num2str(soy) ' sigma0_up = ' num2str(soz) ' ' xresult ' ' yresult ' ' zresult])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % devuelvo los parámetros del ETM y una lista con los elementos que no son outliers
        % devuelvo tambien los pesos actualizados
        [~,oVx]=deleteoutliers(Vx, 0.01);
        [~,oVy]=deleteoutliers(Vy, 0.01);
        [~,oVz]=deleteoutliers(Vz, 0.01);
        outliers = true(size(Vx,1),1);
        outliers(oVx) = 0; outliers(oVy) = 0; outliers(oVz) = 0; 
        
        oindex(i) = num2cell(outliers,1);
        
        weights(i) = num2cell([Px; Py; Pz],[1 2]);
        etm(i,1) = num2cell([Cx'; Cy'; Cz'],[1 2]);
        etm(i,2) = num2cell(Ha,[1 2]);
        
        i = i + 1;
    end
end
