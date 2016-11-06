function [A, Ha] = load_hsf3(stnm, t, lat, Hi)
    % this function builds the design matrix A using the information
    % provided by the H and E files
    
    if nargin > 3
        Ha = Hi;
    else
        Ha = load(['stations/H_' stnm '.txt']);
    end

    [~,its] = min(t);
    [~,ite] = max(t);

    % require at least more than 2 data points to calculate!
    if ite - its > 1
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
            % sparate log decays from jumps
            % Ht has both the decays and jumps to find the offset
            Ht = Ha;

            % heavy side functions
            ht = zeros(size(t,1),size(Ht,1));
            
            % loop through the heaviside functions
            for i = 1:size(Ht,1)
                % check when it starts
                if Ht(i,1) > t(its) & Ht(i,1) < t(ite)
                    % inside time window, put ones
                    ht(:,i) = (t >= Ht(i,1));

                    % if 2 or more jumps happen inside a gap, the (t >= Ht(i,1))
                    % will return identical columns creating a singular matrix.
                    % Verify that the new column is different than the rest
                    if i > 1
                        % if a previous jump exists with less than 3 data points
                        % to constrain it, remove it
                        for j = 1:(i-1)
                            if sum(xor(ht(:,j),ht(:,i))) < 3
                                % identical columns! remove one
                                ht(:,j) = 0;
                                % remove the item from the Ha list
                                Ha(j,:) = 0;
                            end
                        end
                    end
                else
                    % no matching data: remove the item from the Ha list
                    Ha(i,:) = 0;
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
                        % find the index of this jump in Ha
                        k = find(Ha(:,1) == Hl(j,1));
                        
                        if sum(hl(:,j)) == 0 & Ha(k,1) ~= 0
                            % all zeros. Remove from Ha to keep
                            % compatibility with full time serie (ts) vs
                            % actual measured epochs
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