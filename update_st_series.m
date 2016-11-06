function [st_series_o] = update_st_series(st_series,index,weights)
    % update the st_series structure
    
    st_series_o = st_series;

    for i = 1:size(st_series,2)
        
        p = weights{i};
        px = p(1,:);
        py = p(2,:);
        pz = p(3,:);
        
        st_series_o(i).index = index{i};
        
        st_series_o(i).px = px;
        st_series_o(i).py = py;
        st_series_o(i).pz = pz;

        dx = st_series_o(i).x' - mean(st_series_o(i).x);
        dy = st_series_o(i).y' - mean(st_series_o(i).y);
        dz = st_series_o(i).z' - mean(st_series_o(i).z);

        [n,e,u]=ct2lg(dx,dy,dz,st_series_o(i).lat*pi/180,st_series_o(i).lon*pi/180);

        st_series_o(i).n = n';
        st_series_o(i).e = e';
        st_series_o(i).d = u';
    end
end

