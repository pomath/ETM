function st_info = load_st_info(file)
    % read the station.info file
    fileID = fopen(file,'r');

    i = 1;
    tline = fgetl(fileID);
    while ischar(tline)
        if size(tline,2) > 0
            if strcmp(tline(1),' ') && size(tline,2) > 55
                % save the station name, start and end date
                st_info(1,i).stnm = tline(2:5);
                st_info(1,i).ant = tline(171:187);
                st_info(1,i).sdate = str2double(tline(26:29)) + str2double(tline(31:33))/365;
                st_info(1,i).edate = str2double(tline(45:48)) + str2double(tline(50:52))/365;
                i = i + 1;
            end
        end
        tline = fgetl(fileID);
    end

    fclose(fileID);
end