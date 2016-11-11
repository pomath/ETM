clear all
clc

[poly,st_series] = load_polyhedra();

% load the station info
st_info = load_st_info('../tables/station.info');

stab_sites = load_stab_sites('../tables/stab_sites.txt');

adjust_polyhedra(st_series,poly,st_info,stab_sites)

% save('st_series/st_series.mat','st_series');
% save('st_series/etm.mat','etm');
% save('st_series/poly.mat','poly');