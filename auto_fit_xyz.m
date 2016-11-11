function [etm, oindex, weights] = auto_fit_xyz(st_series,st_info,start_at)

    % initialize the cell arrays to save the etm data
    etm = cell(size(st_series,2),2);
    % initialize the cell array for the outlier index information
    oindex = cell(size(st_series,2),1);
    % weights for the GPS data
    weights = cell(size(st_series,2),1);
    
    confidence = 0.05;
    
    % station counter
    if nargin == 3
        i = structfind(st_series,'stnm',start_at);
    else
        i = 1;
    end

    while i <= size(st_series,2)

        stnm = st_series(i).stnm;

        % read in antenna changes from station info and co-seismic jumps
        check_files(stnm, st_info, st_series(i).lat, st_series(i).lon);

        % be aware of any outliers
        index = st_series(i).index;
        
        % get the epochs for the current station
        t = (st_series(i).epochs)';
        
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
        % DDG 10/11/2016: assign an a priori sigma for all observations
        px = ones(size(st_series(i).px));
        py = ones(size(st_series(i).px));
        pz = ones(size(st_series(i).px));
        
        % load the heavy-side functions and periodic terms
        [A,Ha,constrains_h,constrains_s,adjcmp] = load_hsf3(stnm, t,true);
        
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
            factor = 1;
            
            % goodness of fit test
            cst_pass = false;
            % print current component
            fprintf('[%i]',comp);
            
            while ~cst_pass
                % LSQ and goodness of fit adjustment
                % if it doesn't pass, it rewights all the data and takes into account outliers
                % it iterates until it test is ok.
                % should include a max iter
                
                % constrains tells the LSQ where to put the conditions equations to limit the behavior of the log decays and
                % the antenna and co-seismic jumps
                [C, S, So, V, ~, dof, cst_pass] = adjust_lsq(A,P,L,index,constrains_h,constrains_s);

                % verify if the goodness of fit test was passed
                if ~cst_pass
                    if So < 1
                        % weights are too pesimistic, just inform the user
                        fprintf(char(hex2dec('25B2')));
                    else
                        % weights are too optimistic, just inform the user
                        fprintf(char(hex2dec('25BC')));
                    end
                    
                    factor = factor.*So;
                    
                    % determine the outliers
                    [~,outliers]=deleteoutliers(V, confidence);
                    % downweight them (x10)
                    factor_outliers=ones(size(p));
                    factor_outliers(outliers) = 10;
                    % make the new P matrix
                    P = diag(1./((factor.*factor_outliers).^2));
                end
                
                % uncomment this block for debuging purposes
%                 clf
%                 scatter(t,L,30,abs(factor),'fill')
%                 hold on
%                 plot(t,A*C,'r')
%                 colorbar
%                 if ~isempty(outliers)
%                    plot(t(outliers),L(outliers),'xr','MarkerSize',20)
%                 end
            end

            % update the information
            switch comp
                case 1
                    Cx = C; Sx = S; sox = So; Vx = V; Px = p.*factor.*factor_outliers;
                case 2
                    Cy = C; Sy = S; soy = So; Vy = V; Py = p.*factor.*factor_outliers;
                case 3
                    Cz = C; Sz = S; soz = So; Vz = V; Pz = p.*factor.*factor_outliers;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % goodness of fit test
        % this is just to show the value of sigma zero after the test was passed
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        
        % figure out which are outliers and which are not and put zeros in outliers vector
        [~,oVx]=deleteoutliers(Vx, confidence);
        [~,oVy]=deleteoutliers(Vy, confidence);
        [~,oVz]=deleteoutliers(Vz, confidence);
        outliers = true(size(Vx,1),1);
        outliers(oVx) = 0; outliers(oVy) = 0; outliers(oVz) = 0; 
        
        % pack everything in the cell array
        oindex(i) = num2cell(outliers,1);
        
        weights(i) = num2cell([Px; Py; Pz],[1 2]);
        etm(i,1) = num2cell([Cx'; Cy'; Cz'],[1 2]);
        etm(i,2) = num2cell(Ha,[1 2]);
        
        % save_etm_params(adjcmp,min(t),max(t),[Cx Cy Cz],Ha,stnm)
        i = i + 1;
    end
end
