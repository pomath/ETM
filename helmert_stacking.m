function [st_series_o, etm, poly] = helmert_stacking(st_series, poly, st_info, etm, iter, stab_sites)

    X = poly.x;
    Y = poly.y;
    Z = poly.z;
    t = poly.epochs;
    
    % numero de estaciones
    n = size(X,2);
    % numero de epocas
    e = size(poly.epochs,1);
    
    % aqui esta la iteracion principal de la rotacion y traslacion de poliedros
    % se utilizan todas las estaciones
    
    rot_hist = cell(e,iter);
    
    for j = 1:iter

        % vector of epochs to delete
        del_epochs = [];
        
        % calculate the expected coordinates from the ETM
        [Xe, Ye, Ze] = get_coords_etm(st_series,etm,t);
        
        % identify outliers and stabilization station
        out = get_outliers_stab(st_series,poly,stab_sites);
    
        % matrix to put the results
        Lr = nan(e,3*n);
        
        for i = 1:e
            
            % matriz de diseño para cada componente
            Ax = [zeros(n,1) -Z(i,:)' Y(i,:)' X(i,:)' repmat([1 0 0],n,1)];
            Ay = [Z(i,:)' zeros(n,1) -X(i,:)' Y(i,:)' repmat([0 1 0],n,1)];
            Az = [-Y(i,:)' X(i,:)' zeros(n,1) Z(i,:)' repmat([0 0 1],n,1)];
            %Ax = [zeros(n,1) -Z(i,:)' Y(i,:)' repmat([1 0 0],n,1)];
            %Ay = [Z(i,:)' zeros(n,1) -X(i,:)' repmat([0 1 0],n,1)];
            %Az = [-Y(i,:)' X(i,:)' zeros(n,1) repmat([0 0 1],n,1)];
            
            % vector de observaciones
            L = [Xe(:,i) - X(i,:)'; Ye(:,i) - Y(i,:)'; Ze(:,i) - Z(i,:)'];

            % junto las tres componentes
            A = [Ax(:,:); Ay(:,:); Az(:,:)];

            % remove rows with zeros (those that don't have observations)
            % la estructura de poliedros exige que todos los vectores tengan el
            % mismo tamaño. Algunas epocas no tienen datos, porque no se pudo
            % resolver en GAMIT o porque el dato fue descartado por la primera
            % pasada de auto_fit_xyz. Detecto estos faltantes porque los
            % valores estan en NaN.
            A([out(:,i); out(:,i); out(:,i)] == 0,:) = [];
            L([out(:,i); out(:,i); out(:,i)] == 0,:) = [];

            % build a unitary weights matrix
            p = ones(size(L));
            P = diag(p);
            
            % si el sistema esta sobredeterminado, resuelvo
            if size(A,1) > size(A,2)

                [x, ~,So,~,~,~,cst_pass,~,index] = adjust_lsq(A,P,L,[],[]);
            
                if cst_pass
                    cst_pass = 'PASS';
                else
                    cst_pass = 'NOT PASS';
                end
                
                if sqrt(x(end-2).^2 + x(end-1).^2 + x(end).^2) > 0.003;
                    fprintf('\n %4.3f (%3.2f: %s) trans (mm): [%+4.3f %+4.3f %+4.3f] >> outliers detected: ',t(i), So, cst_pass,x(end-2).*1000,x(end-1).*1000,x(end).*1000)
                    fprintf('%d ', find(index == 0))    
                else
                    fprintf('\n %4.3f (%3.2f: %s); max disp: %3.2f mm',t(i), So, cst_pass, max(abs(A*x)).*1000)
                end
                fprintf('\n')
                
                % make A again to include the all stations back (
                A = [Ax; Ay; Az];

                % apply the rotation and translation and scale factor
                Lr(i,:) = (A*x)' + [X(i,:) Y(i,:) Z(i,:)];

                % the new X Y Z coordinates to perform the next adjustment are the rotated ones
                % cuando terminen todas las epocas, se vuelve a emplear el ETM
                % para ver cuanto hay que rotar. Se guarda en X Y Z las
                % coordenadas de los poliedros ya rotados y trasladados.
                X(i,:) = (Ax*x)' + X(i,:);
                Y(i,:) = (Ay*x)' + Y(i,:);
                Z(i,:) = (Az*x)' + Z(i,:);

                rot_hist(i,j) = num2cell([t(i); x],1);
            else
                % insufficient observation equations!
                % remove points from time series
                % save a vector of the epochs to delete
                del_epochs = [del_epochs; poly.epochs(i)];
            end
        end

        % rebuild the st_series from the Lr matrix
        st_series = load_poly2series(st_series,Lr,del_epochs);

        % refit without common mode
        % DDG: el ajuste y obtención de los parametros para la remocion de modo
        % comun no esta implementada todavia. Esto es porque ahora cada
        % estacion tiene un juego de parametros diferente cuando las
        % componentes de marea no son iguales
        if j < iter
            [etm, index, weights] = auto_fit_xyz(st_series,st_info,stab_sites);
            st_series = update_st_series(st_series,index,weights,stab_sites);
        end
        
        disp(['Run ' num2str(j) ' of ' num2str(iter)])
    end

    % save the roto-translation history
    save('../ETM_files/st_series/rot_hist.mat','rot_hist');
    
    % now that we have the st_series aligned, we can adjust an ETM to all
    % the stations (not only the stab_sites)
    [etm, index, weights] = auto_fit_xyz(st_series,st_info);
    st_series = update_st_series(st_series,index,weights);
    
    % rebuild the polyhedron struct
    poly = create_poly_struct(st_series);
    
    st_series_o = st_series;
end

