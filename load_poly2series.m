function [st_series_o]= load_poly2series(st_series,Lr,del_epochs)

    stnc = size(st_series,2);
    st_series_o(stnc) = struct();
    for i = 1:size(st_series,2)
        
        % copy the data that is identical
        for fn = fieldnames(st_series)'
            st_series_o(i).(fn{1}) = st_series(1,i).(fn{1});
        end
        st_series_o(i).x = Lr(:,i)';
        st_series_o(i).y = Lr(:,stnc+i)';
        st_series_o(i).z = Lr(:,stnc*2+i)';
        
        % delete coordinates with NaNs (probably bad data or could not
        % solve system of equations)
        d = sqrt(st_series_o(i).x.^2 + st_series_o(i).y.^2 + st_series_o(i).z.^2);
        
        %st_series_o(1,i).x(d < 6300e3) = [];
        st_series_o(i).x(isnan(st_series_o(i).x)) = [];
        %st_series_o(1,i).y(d < 6300e3) = [];
        st_series_o(i).y(isnan(st_series_o(i).y)) = [];
        %st_series_o(1,i).z(d < 6300e3) = [];
        st_series_o(i).z(isnan(st_series_o(i).z)) = [];
        
        st_series_o(i).n(ismember(st_series_o(i).epochs,del_epochs)) = [];
        st_series_o(i).e(ismember(st_series_o(i).epochs,del_epochs)) = [];
        st_series_o(i).d(ismember(st_series_o(i).epochs,del_epochs)) = [];
        st_series_o(i).pn(ismember(st_series_o(i).epochs,del_epochs)) = [];
        st_series_o(i).pe(ismember(st_series_o(i).epochs,del_epochs)) = [];
        st_series_o(i).pd(ismember(st_series_o(i).epochs,del_epochs)) = [];
        st_series_o(i).px(ismember(st_series_o(i).epochs,del_epochs)) = [];
        st_series_o(i).py(ismember(st_series_o(i).epochs,del_epochs)) = [];
        st_series_o(i).pz(ismember(st_series_o(i).epochs,del_epochs)) = [];
        st_series_o(i).index(ismember(st_series_o(i).epochs,del_epochs)) = [];
        
        st_series_o(i).epochs(ismember(st_series_o(i).epochs,del_epochs)) = [];
    end
end
