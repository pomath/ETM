function [A, Ha, constrains_h, constrains_s, adjcmp] = load_hsf3(stnm, t, inversion, Hi)
    % this function builds the design matrix A using the information
    % provided by the H and E files
    
    if nargin > 3
        Ha = Hi;
    else
        Ha = load(['../ETM_files/stations/H_' stnm '.txt']);
    end

    [~,its] = min(t);
    [~,ite] = max(t);

    % require at least more than 2 data points to calculate!
    if ite - its > 0
        % determine the start of this sections of time series
        Tr = t(its);

        % constant offset variable
        c = ones(size(t,1),1);
        
        % velocity variable
        v = (t - Tr);
        
        s = [sin(2*pi*t) sin(4*pi*t) cos(2*pi*t) cos(4*pi*t)];
        
        a = [];
        
        if ~isempty(Ha)
            % reorder Ha
            Ha = sortrows(Ha,1);
            
            % remove jumps (no log decay) before and after the time
            % limits. Those don't influence the design matrix
            Ha(Ha(:,2) == 0 & (Ha(:,1) <= t(its) | Ha(:,1) >= t(ite)),:) = [];
            
            % remove any log decays that happened after the t(ite)
            Ha(Ha(:,2) ~= 0 & Ha(:,1) >= t(ite),:) = [];
            
            % log decays < t(its) are more complicated. There are a few
            % possibilities:
            % 1) The jump happened before the data started, but the decay
            % will continue. The decay needs to be added but not the jump
            % 2) Idem before but the decay is stopped by another decay that
            % happened < t(its) => remove all together
            
            % find the min(t(its) - t_decay) but limit the search to Ha <= t(its) and ...
            min_Hl = t(its) - min(t(its) - Ha(Ha(:,2) ~= 0 & Ha(:,1) <= t(its)));
            if ~isempty(min_Hl)
                % remove everything < than t(its) - min(t(its) - t_decay)
                Ha(Ha(:,2) ~= 0 & Ha(:,1) < min_Hl,:) = [];
            end
            
            % heavy side functions
            ht = zeros(size(t,1),size(Ha,1));
            
            % loop through the heaviside functions
            for i = 1:size(Ha,1)
                % check when it starts
                % the jump has to be within t(its) and t(ite) to get a jump
                if Ha(i,1) > t(its) & Ha(i,1) < t(ite)
                    % inside time window, put ones
                    ht(:,i) = (t >= Ha(i,1));

                    % if 2 or more jumps happen inside a gap, the (t >= Ht(i,1))
                    % will return identical columns creating a singular matrix.
                    % Verify that the new column is different than the rest
                    if i > 1
                        % if a previous jump exists with less than 2 data points
                        % to constrain it, remove it
                        for j = 1:(i-1)
                            % compare col j with col i
                            if sum(xor(ht(:,j),ht(:,i))) < 2
                                % identical columns! remove the previous one (only
                                % different by 2 elements)
                                ht(:,j) = 0;
                                % remove the item from the Ha list
                                Ha(j,:) = 0;
                            end
                        end
                    end
                end
            end
            
            Hl = Ha(Ha(:,2) ~= 0,:);
            hl = zeros(size(t,1),size(Hl,1));
            
            % loop through the heaviside functions (log decays)
            for i = 1:size(Hl,1)
                hl(:,i) = (t > Hl(i,1)).*log10(1+(t-Hl(i,1))./Hl(i,2));

                % check for previous log decays. Stop them if necessary
                if i > 1
                    % if a previous decay exists, stop it
                    for j = 1:(i-1)
                        hl(t >= Hl(i,1),j) = 0;
                        if sum(hl(:,j)) == 0
                            % if the col j is all zeros, remove this log
                            % decay (probably no data to constrain it)
                            % remove also from Ha. This is because the
                            % inversion will not constrain this decay and
                            % when we call load_hsf3 with a full time
                            % vector (named ts) it will try to put it and
                            % we'll end with one more column than
                            % parameters in the inversion.
                            
                            % find the index of this jump in Ha
                            k = Ha(:,1) == Hl(j,1);
                            Ha(k,:) = 0;
                        end
                    end
                end
            end       

            % verify that all cols are not zero
            ht(:,~any(ht,1)) = [];
            hl(:,~any(hl,1)) = [];
            
            Ha(~any(Ha,2),:) = [];
        else
            ht = [];
            hl = [];
        end
        
        A = [c v a ht hl s];
        if inversion == true
            % add some constrains to the A matrix to stabilize the
            % inversion of the log decays and jumps
            % constrains = [zeros(size([ht hl],2),size([c v a],2)) diag(ones(size([ht hl],2),1)) zeros(size([ht hl],2),size(s,2))];
            constrains_h = [zeros(size([ht hl],2),size([c v a],2)) diag(ones(size([ht hl],2),1)) zeros(size([ht hl],2),size(s,2))];
            constrains_s = [zeros(size(s,2),size([c v a ht hl],2)) diag(ones(size(s,2),1))];
        end
        
        adjcmp = [size(c,2) size(v,2) size(a,2) size(ht,2) size(hl,2) size(s,2)];
    end
end

%         m(1) = 2*pi/((12.4206/24)/365);
%         m(2) = 2*pi/((12.0000/24)/365);
%         m(3) = 2*pi/((12.6584/24)/365);
%         m(4) = 2*pi/((11.9673/24)/365);
%         m(5) = 2*pi/((23.9344/24)/365);
%         m(6) = 2*pi/((25.8194/24)/365);
%         m(7) = 2*pi/((24.0659/24)/365);
%         m(8) = 2*pi/((5.9836/24)/365);
%         m(9) = 2*pi/(13.661/365);
%         m(10) = 2*pi/(27.555/365);
% 
%         s = [];
%         for i = 1:length(m)
%             s = [s sin(m(i)*t) cos(m(i)*t)];
%         end