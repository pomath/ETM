function [st_series_o, etm] = helmert_tides(st_series, poly, st_info, etm, iter)

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
    for j = 1:iter

        % matrix to put the results
        Lr = [];
        % vector of epochs to delete
        del_epochs = [];
        
        % calculate the expected coordinates from the ETM
        [Xe, Ye, Ze] = get_coords_etm(st_series,etm,t);
            
        for i = 1:e
            
            % identify outlier stations for this epoch
            out = get_outliers(st_series,t(i));
            
            % matriz de diseño para cada componente
            Ax = [zeros(n,1) Z(i,:)' -Y(i,:)' repmat([1 0 0],n,1)];
            Ay = [-Z(i,:)' zeros(n,1) X(i,:)' repmat([0 1 0],n,1)];
            Az = [Y(i,:)' -X(i,:)' zeros(n,1) repmat([0 0 1],n,1)];

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
            A([out; out; out] == 0,:) = [];
            L([out; out; out] == 0,:) = [];

            % si el sistema esta sobredeterminado, resuelvo
            if size(A,1) > size(A,2)

                % sin pesos por el momento
                x = (A'*A)\(A'*L);

                % make A again to include the all stations back (
                A = [Ax; Ay; Az];

                % apply the rotation and translation
                Lr = [Lr; (A*x)' + [X(i,:) Y(i,:) Z(i,:)]];

                % the new X Y Z coordinates to perform the next adjustment are the rotated ones
                % cuando terminen todas las epocas, se vuelve a emplear el ETM
                % para ver cuanto hay que rotar. Se guarda en X Y Z las
                % coordenadas de los poliedros ya rotados y trasladados.
                X(i,:) = (Ax*x)' + X(i,:);
                Y(i,:) = (Ay*x)' + Y(i,:);
                Z(i,:) = (Az*x)' + Z(i,:);

                rot_hist(i,j) = num2cell(x,1);
            else
                % insufficient observation equations!
                % remove points from time series
                Lr = [Lr; nan(1,size(X,2).*3)];
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
        [etm, index, weights] = auto_fit_xyz(st_series,st_info);
        
        st_series = update_st_series(st_series,index,weights);

        disp(['Run ' num2str(j) ' of ' num2str(iter)])
    end

    st_series_o = st_series;
end

