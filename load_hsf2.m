function [A, Ha, Ec, fits_out, Hsf, ldc_out, dt_out] = load_hsf2(stnm, t, overide, from_model)
    % this function builds the design matrix A using the information
    % provided by the H and E files
    
    Ha = load(['stations/H_' stnm '.txt']);
    Ec = load(['stations/E_' stnm '.txt']);
    
    % everything is masked by the E file cut
    % define as many linear fits as lines the E file has.
    % take into account the size of t, so we dont't put any unnecessary
    % fits (where there is no data)
    
    fits = size(Ec,1);
    fits_index = zeros(fits,1);
    fits_out = 0;
    c = [];
    v = [];
    s = [];
    a = [];
    h = [];
    dt_out = [];
    Ht= [];
    
    for f = 1:fits
        % determine t start t end
        % choose between the established min or the min vector value
        vts = Ec(f,1);
        if f == fits && max(t) > Ec(f,2) && overide
            vte = max(t);
        else
            vte = Ec(f,2);
        end
        
        % check if there is data inside the time window
        if min(t) <= vte && max(t) >= vts
            
            its = find(t >  vts, 1, 'first');
            ite = find(t <= vte, 1, 'last');

            % require at least more than 2 data points to calculate!
            if ite - its > 1
                fits_index(f) = true;
                
                fits_out = fits_out + 1;
                % determine the start of this sections of time series
                Tr = t(its);

                % constant offset variable
                c = [c (t>=t(its) & t<=t(ite))];
                % velocity variable
                if from_model == 1
                    v = [];
                else
                    v = [v (t - Tr)];
                    % make zero everything outside the time window
                    v(t<t(its),end) = 0;
                    v(t>t(ite),end) = 0;
                end
                
                s = [s sin(2*pi*t) cos(2*pi*t) sin(4*pi*t) cos(4*pi*t)];
                s(t<t(its),end-3:end) = 0;
                s(t>t(ite),end-3:end) = 0;
                
                a = [a (t - Tr).^2];
                a(t<t(its),end) = 0;
                a(t>t(ite),end) = 0;
                
                % heavy side functions
                ldc = 1; % log decay counts
                ldc_out(fits_out) = 0;

                hti = 1; % Ht index
                one_trans = 0;
                for j = 1:size(Ha,1)
                    % if the jump in the H file is within this time period, add it
                    if Ha(j,1) >= vts & Ha(j,1) < vte & sum(t >= Ha(j,1)) > 0 & ...
                            (sum(t >= Ha(j,1)) ~= size(t,1) | one_trans == 0)

                        Ht(:,hti) = (t >= Ha(j,1));
                        
                        % if there is a transient for the whole time
                        % series, do not allow more then one (that fills
                        % the whole time series)
                        if sum(t >= Ha(j,1)) == size(t,1) 
                            one_trans = 1;
                        end
                        
                        % check if Ht resulted in all zeros. This happens when
                        % the time series stop before an antenna change was
                        % reported. Eg: CFAG antenna changed in 2012.5 but time
                        % series stop in 2011, so t >= Ha(j,1) returns all 0
                        % ALSO, the step function for antenna changes only
                        % makes sense if there is data before and after! if
                        % there is no data before and after, then there is no
                        % point in having the step function, which will make
                        % the A matrix singular (because of the columns of ones
                        % to find the offset)

                        % make zero everything outside the time window
                        if sum(Ht(:,hti)) == 0 || (Ha(j,2) == 0 && (Ha(j,1) <= t(its) || Ha(j,1) >= t(ite)))
                            Ht(:,hti) = [];
                        else
                            Ht(t<t(its),hti) = 0;
                            Ht(t>t(ite),hti) = 0;
                            hti = hti + 1;
                        end

                        % if column 2 ~= 0, then we need to apply a log after the Ht
                        if Ha(j,2) ~= 0
                            % save the indecies where a log decay is needed
                            Hl(:,ldc) = (t > Ha(j,1));

                            % make zero everything outside the time window
                            Hl(t<vts,ldc) = 0;
                            Hl(t>vte,ldc) = 0;

                            % time of the earthquake
                            dT(1,ldc) = Ha(j,1);
                            % relaxation time
                            T(1,ldc) = Ha(j,2);
                            dt_out = T;

                            % make the jump on Ht also zero, since this is an
                            % earthquake jump and the log decay forces to make a
                            % new fit => the jump is taken care by the linear
                            % component
                            % NO LONGER TRUE AFTER SERGIO'S DECISION OF
                            % USING THE FULL BEVIS AND BROWN DEFINITION
                            % salvo que no haya datos anteriores al sismo
                            % entonces, hay 2 columnas de unos (1), una por
                            % la componente lineal (offset) y otra por la
                            % componente del salto.
                            if Ha(j,1) <= t(its)
                                % como el comienzo de la serie de tiempo es
                                % posterior a la ubicacion del salto, hay
                                % que sacar la columna de todos unos.
                                Ht(:,hti-1) = [];
                                hti = hti - 1;
                            end
                            
                            % if there is more than one decay, stop the previous fit
                            % to start with the new one.
                            if ldc > 1
                                Hl(:,ldc-1) = Hl(:,ldc-1) & ~Hl(:,ldc);
                            end
                            ldc = ldc+1;
                            ldc_out(fits_out) = ldc_out(fits_out) + 1;
                        end
                    end
                end

                % build the A columns for the log decay
                if ldc > 1
                    h = [h Hl.*log10(1+(repmat(t,1,size(dT,2))-repmat(dT,[size(t),1]))./repmat(T,[size(t),1]))];
                end
            end
        end
    end
    
	A = [c v Ht h s];
    Hsf = size(Ht,2);
    Ec = Ec(logical(fits_index),:);
end