clear all
clc

[poly,st_series] = load_polyhedra();

% load the station info
% el station info se usa para agregar saltos por cambios de equipo, antena,
% etc.
st_info = load_st_info('tables/station.info_demian');

stab_sites = load_stab_sites('tables/stab_sites_demian.txt');

