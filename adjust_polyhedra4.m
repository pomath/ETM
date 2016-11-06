% DDG: 12/06/2016
% script para ajustar la marea terrestre utilizando el ETM y series de
% tiempo GPS. Algunos comentarios están en inglés porque los archivos
% provienen de scripts que vengo escribiendo para la doctorado.

%{

Run make_data() then this.

%}

function adjust_polyhedra4();

% load the time series
[poly,st_series] = load_polyhedra();

% load the station info
% el station info se usa para agregar saltos por cambios de equipo, antena,
% etc.
st_info = load_st_info('/home/pmatheny/scot/tables/station.info');

% fit etms and filter large outliers
% primera corrida de auto_fit_xyz busca datos que estén por afuera de 3
% sigma y los descarta para no deformar mucho la red durante el ajuste
[etm, index, weights] = auto_fit_xyz(st_series,st_info);

% actualizo las series de tiempo sacando las colas estadisticas
st_series = update_st_series(st_series,index,weights);

% [st_series, etm] = helmert_tides(st_series, poly, st_info, etm, 4);

plot_etms(st_series, etm)

save('st_series/st_series.mat','st_series');
save('st_series/etm.mat','etm');
save('st_series/poly.mat','poly');
end
