function [ref_sites,ref_sites_i, ref_sites_v] = load_ref_vel(st_series,etm,sinex)
    fid = fopen(sinex,'r');
    
    [sol,~,stnvel_idx] = sinex_read_solutionestimate(fid);
    
    ref_sites = unique(sol{1,3}(stnvel_idx,:));
    
    ref_sites_i = zeros(size(ref_sites,1),1);
    ref_sites_v = zeros(size(ref_sites,1),3);
    for i = 1:length(ref_sites)
        index = find(strcmp({st_series.stnm}, ref_sites{i})==1);
        if ~isempty(index)
            % check if this can be a stabilization site (need to be able to
            % have an ETM and > 50% of the data for its timespan)
            if check_stn_data(st_series(index).epochs,0.5) & ~isempty(etm{index,1})
                % if station has only a post-seismic estimate of the
                % velocity, then disregard. This is because the linear
                % portion will be unstable and will not reflect the real
                % inter-seismic behavior
                Ha = etm{index,2};
                if isempty(Ha(Ha(:,2) ~=0,1))
                    minpost=99999;
                else
                    minpost=Ha(Ha(:,2) ~=0,1);
                end
                if min(st_series(index).epochs) < minpost
                    % find this site in the solution vector
                    soli = strcmp(sol{1,3}, ref_sites{i})==1 & stnvel_idx == 1;
                    if sum(soli) == 3
                        % more than one solution for this site (???)
                        % do not use
                        ref_sites_v(i,:) = sol{1,11}(soli)';
                        ref_sites_i(i) = index;
                    end
                end
            end
        end
    end
    ref = [ref_sites_i ref_sites_v];
    ref = sortrows(ref,1);
    ref(ref(:,1) == 0,:) = [];
    ref_sites_i = ref(:,1);
    ref_sites_v = ref(:,2:4);
    
end