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
        st_series_o(i).x(isnan(st_series_o(i).x)) = [];
        st_series_o(i).y(isnan(st_series_o(i).y)) = [];
        st_series_o(i).z(isnan(st_series_o(i).z)) = [];
        
        st_series_o(i).px(ismember(st_series_o(i).epochs,del_epochs)) = [];
        st_series_o(i).py(ismember(st_series_o(i).epochs,del_epochs)) = [];
        st_series_o(i).pz(ismember(st_series_o(i).epochs,del_epochs)) = [];
        st_series_o(i).index(ismember(st_series_o(i).epochs,del_epochs)) = [];
        
        st_series_o(i).epochs(ismember(st_series_o(i).epochs,del_epochs)) = [];
    end
end
