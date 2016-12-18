function [] = check_files(stnm,st_info,lat,lon)
    % check the existance of the H and E file and create them if not
    
    % if the file does not exist,
    % create the file with the Maule event
    th = [];
    if ~exist(['../ETM_files/stations/H_' stnm '.txt'], 'file')
        % find the antenna changes in the station info file
        %idx = structfind(st_info,'stnm',stnm);
        idx = find(strcmp({st_info.stnm}, stnm)==1);
        atx_H = [];
        if ~isempty(idx)
            atx_sdate = [st_info(idx).sdate];
            atx_edate = [st_info(idx).edate];
            % save the antennas in a cell array
            antennas = {st_info(idx).ant};
            
            % resort the rows just in case
            atx_date = sortrows([atx_sdate' atx_edate' (1:size(atx_edate,2))'],[1 2]);
            if size(atx_date,1) > 1
                for h = 2:size(atx_date,1)
                    atx_H = [atx_H; atx_date(h,1)];
                end
            end
            % keep unique dates
            atx_H = [unique(atx_H) zeros(size(unique(atx_H),1),1)];
        end
        
        % open the file of seismic events
        fileID = fopen('../ETM_files/seismic/events.csv','r');
        % skip the header
        fgetl(fileID);
        events = textscan(fileID, '%s%f%f%f%s%s%f%f%f%s%s%s%s%q%s%s%s%s%s%s%s%s','delimiter',',');
        
        edist = distance(lat,lon,events{2},events{3})*pi/180*6371;
        
        
        mw = events{5};
        tS = sprintf('%s ', mw{:});
        mw = sscanf(tS, '%f');
        
        % Mike's expression
        S = -0.8717*(log10(edist)-2.25) + 0.4901*(mw-6.6928);
        
        % magn = (mw-6.5)*700;
        % FROM NEVADA'S WEBSITE:
        % magn = 10.^(0.5.*mw - 0.8);
        
        edepth = events{4};
        
        edate = datevec(datenum(strrep(strrep(events{1},'T',' '),'Z', ''),'yyyy-mm-dd HH:MM:ss.FFF'));
        
        %for i = 1:size(edate,1)
            %[doy(i,1),fyear(i,1)] = date2doy(datenum(sprintf('%d/%d/%d %d:00',edate(i,3),edate(i,2),edate(i,1),edate(i,4))));
        %end        
        fyear = date2fyear(edate(:,1),edate(:,2),edate(:,3),edate(:,4)+edate(:,5)./60+edate(:,6)./3600);
        
        evinfo = [fyear edepth edist mw S];
        evinfo = sortrows(evinfo,1);
        % remove events with depth > 50 km
        % evinfo(evinfo(:,2) >= 50,:) = [];
        
        % keep relevant events to this station
        events = evinfo(evinfo(:,5) > 0,:);
        
        if ~isempty(events)
            finished = false;
            if size(events,1) > 1
                while ~finished
                    % check for same day or consecutive events and combine
                    repevents = diff(events(:,1));

                    for i = 1:size(repevents,1)
                        if repevents(i) <= 21/(365+leapyear(floor(events(i,1))))
                            % these two events happened very close to each other in
                            % space and time. Keep the larger event
                            if events(i,5) < events(i+1,5)
                                events(i,:) = [];
                            else
                                events(i+1,:) = [];
                            end
                            break;
                        end
                    end
                    if i == size(repevents,1)
                        finished = true;
                    end
                end
            end

            th = [atx_H; events(:,1) repmat(0.5,size(events,1),1)];
        else
            th = atx_H;
        end
        % Take into account that GPS days starts at DOY = 1
        % add the antena jumps to the co-seismic jump
        % verify the distance to the earthquake and put the jump in of
        % applies
        %if distance(lat,lon,-36.14,-72.93)*pi/180*6371 <= 1500
        %    th = [atx_H; 2010 + 57/365 + 6/(24*365) + 34/(1440*8736) 1];
        %end
        save(['../ETM_files/stations/H_' stnm '.txt'], 'th','-ascii'); 
    end

end