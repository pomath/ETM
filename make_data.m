function make_data();
%{
Run the shell script to generate the time_series.txt file from the globk
output then run this.

%}

fileID = fopen('/home/pmatheny/scripts/time2.txt','r');

A = textscan(fileID,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f');

fclose(fileID);

% slip each station into its own vector
stnv = char(A{1,1});
data = [A{1,2} A{1,3} A{1,4} A{1,5} A{1,6} A{1,7} A{1,8} A{1,9} A{1,10} A{1,11} A{1,12} A{1,13} A{1,14}];
% find all unique GPS station
stns = unique(stnv,'rows');
% build different variables with the data
for i = 1:size(stns,1)
    index = strmatch(stns(i,:),stnv);
    
    if isstrprop(stns(i,1),'digit')
        stname = ['i' stns(i,2:4)];
        eval([stname '= data(index,:);']);
    else
        stname = stns(i,1:4);
        if ~isempty(strfind(stname,'@')); stname = strrep(stname, '@', 'a'); end;
        if ~isempty(strfind(stname,'!')); stname = strrep(stname, '!', 'a'); end;
        eval([stname '= data(index,:);']);
    end
    temp = data(index,:);
    % save the data in text files
    fileID = fopen(['series_new/' stname '.txt'],'w');
    for j = 1:size(temp,1)
        %if sum(temp(j,11:13)) < 0.3
        if jd2yr(mjd2jd(temp(j,1))) > 2008
            % do not allow solutions with errors > 0.3
            fprintf(fileID,'%s %f %+f %+f %+f %+f %+f %+f %+f %+f %+f %+f %+f %+f\n',stname,jd2yr(mjd2jd(temp(j,1))),temp(j,2),temp(j,3),temp(j,4),temp(j,5),temp(j,6),temp(j,7),temp(j,8),temp(j,9),temp(j,10),temp(j,11),temp(j,12),temp(j,13));
        end
    end
    fclose(fileID);
end
end

