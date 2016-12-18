function [st_series_o] = update_st_series(st_series,index,weights,start_at)
    % update the st_series structure
    
    st_series_o = st_series;

    % station counter
    if nargin == 4
        if size(start_at,1) == 1
            j = (find(strcmp({st_series.stnm}, start_at)==1):length(st_series))';
        else
            j = sort(start_at);
        end
    else
        j = (1:length(st_series))';
    end
    
    for i = j'
        
        if ~isempty(weights{i})
            p = weights{i};
            px = p(1,:);
            py = p(2,:);
            pz = p(3,:);

            st_series_o(i).index = index{i};

            st_series_o(i).px = px;
            st_series_o(i).py = py;
            st_series_o(i).pz = pz;
        end
        
        dx = st_series_o(i).x' - mean(st_series_o(i).x);
        dy = st_series_o(i).y' - mean(st_series_o(i).y);
        dz = st_series_o(i).z' - mean(st_series_o(i).z);

        [n,e,u]=ct2lg(dx,dy,dz,st_series_o(i).lat*pi/180,st_series_o(i).lon*pi/180);

        st_series_o(i).n = n';
        st_series_o(i).e = e';
        st_series_o(i).d = u';

    end
end

