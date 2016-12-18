function [Xe,Ye,Ze] = get_coords_etm(st_series,etm,ts)
% get_coords_etm calculates the ETM predicted coordinates (X,Y,Z) for stations in st_series
% parámetros estimados del ETM
%    
    Xe = nan(size(st_series,2),size(ts,1));
    Ye = nan(size(st_series,2),size(ts,1));
    Ze = nan(size(st_series,2),size(ts,1));

    for i = 1:size(st_series,2)
        stnm = st_series(i).stnm;
        Ha = etm{i,2};
        
        % do not load the Ha from file; use the Ha from the real data
        % adjustment
        % do not send ts in full. Will cause problems with heavy side
        % functions from log decay that occured before there was data in
        % the time series only calculate A with data within limits and fill
        % the rest with NaNs
        t = st_series(i).epochs';
%     
%         [A,~] = load_hsf3(stnm, t, false,Ha);
% 
%         C = etm{i,1};
%         Cx = C(1,:)';
%         Cy = C(2,:)';
%         Cz = C(3,:)';
%         
%         if size(A,2) ~=size(Cx,1)
%             disp(['Warning: ' stnm ' data and A matrix size do not agree'])
%         end
%         
        %tindex = ismember(ts,t);
        
        %Xe(i,tindex) = A*Cx;
        %Ye(i,tindex) = A*Cy;
        %Ze(i,tindex) = A*Cz;
        
        tst = -1;
        if ~isempty(Ha)
            for j = 1:size(Ha,1)

                tindex = ismember(ts,t(t>=tst & t<Ha(j,1)));

                if sum(tindex) > 60
                    w = 60;
                else
                    w = ceil(sum(tindex)/4);
                end
                if w > 1
                    [y,~,~]=ssa(st_series(i).x(t>=tst & t<Ha(j,1))',w,1);
                    Xe(i,tindex) = y;
                    [y,~,~]=ssa(st_series(i).y(t>=tst & t<Ha(j,1))',w,1);
                    Ye(i,tindex) = y;
                    [y,~,~]=ssa(st_series(i).z(t>=tst & t<Ha(j,1))',w,1);
                    Ze(i,tindex) = y;
                else
                    Xe(i,tindex) = st_series(i).x(t>=tst & t<Ha(j,1))';
                    Ye(i,tindex) = st_series(i).y(t>=tst & t<Ha(j,1))';
                    Ze(i,tindex) = st_series(i).z(t>=tst & t<Ha(j,1))';
                end
                tst = Ha(j,1);
            end
            disp(stnm)
            tindex = ismember(ts,t(t>=Ha(end,1)));

            if sum(tindex) > 60
                w = 60;
            else
                w = ceil(sum(tindex)/4);
            end
            if w > 1
                [y,~,~]=ssa(st_series(i).x(t>=Ha(end,1))',w,1);
                Xe(i,tindex) = y;
                [y,~,~]=ssa(st_series(i).y(t>=Ha(end,1))',w,1);
                Ye(i,tindex) = y;
                [y,~,~]=ssa(st_series(i).z(t>=Ha(end,1))',w,1);
                Ze(i,tindex) = y;
            else
                Xe(i,tindex) = st_series(i).x(t>=Ha(end,1))';
                Ye(i,tindex) = st_series(i).y(t>=Ha(end,1))';
                Ze(i,tindex) = st_series(i).z(t>=Ha(end,1))';
            end
        else
            if sum(length(st_series(i).x)) > 60
                w = 60;
            else
                w = ceil(length(st_series(i).x)/4);
            end
            
            tindex = ismember(ts,t);
            if w > 1
                [y,~,~]=ssa(st_series(i).x',w,1);
                Xe(i,tindex) = y;
                [y,~,~]=ssa(st_series(i).y',w,1);
                Ye(i,tindex) = y;
                [y,~,~]=ssa(st_series(i).z',w,1);
                Ze(i,tindex) = y;
            else
                Xe(i,tindex) = st_series(i).x';
                Ye(i,tindex) = st_series(i).y';
                Ze(i,tindex) = st_series(i).z';
            end
        end
        if strcmp(stnm,'CONZ') == 1
           disp('sant') 
        end
    end
end
