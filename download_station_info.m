% SCRIPT TO download the station info files from UNAVCO
% awk 'FNR==NR {a[$0]++; next} !a[$0]' station.info.mit station.info.ign > station.info.ign.only
% 
% stab_sites = load_stab_sites('../tables/stab_sites.txt');
% 
% for i = 1:length(stab_sites)
%     unix(['curl http://www.unavco.org/data/web-services/gsacws/gsacapi/site/search/sites?output=site.station.info\&site.code=' stab_sites{i} ' | grep -v "9999 999 00 00 00  9999 999 00 00 00" > ../tables/stn_info_new/' stab_sites{i} '.TXT'])
%     unix(['cat ../tables/stn_info_new/' stab_sites{i} '.TXT >> ../tables/stn_info_new/station.info'])
% end

files=dir('../tables/stn_info_new/*.log');

for i = 1:length(files)
    unix(['sh_upd_stnfo -ref ../tables/station.info -i ../tables/stn_info_new/' files(i).name '.log'])
end