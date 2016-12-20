clear all
clc

[poly,st_series] = load_polyhedra('../ETM_files/series');

% load the station info
st_info = load_st_info('../tables/station.info.new');

[stab_sites,stab_sites_i] = load_stab_sites(st_series,'../tables/stab_sites.txt');

adjust_polyhedra(st_series,poly,st_info,stab_sites_i)

% plot the stabilization sites
load('../ETM_files/st_series/stab_sites_i.mat');
load('../ETM_files/st_series/st_series.mat');
load('../ETM_files/st_series/etm.mat');

plot_etms(st_series, etm, 'etm')

% plot the stab_sites
clf
subplot(4,2,[1 3 5 7])
m_proj('mercator')
m_coast();
hold on
lon=[st_series(stab_sites_i).lon]';
lon(lon > 180) = lon(lon > 180) - 360;
lat=[st_series(stab_sites_i).lat]';
stnm=char({st_series(stab_sites_i).stnm});
m_plot(lon,lat,'or');
m_text(lon,lat,stnm,'FontSize',6);
axis equal
axis tight

% plot the rotation-translation history
load('../ETM_files/st_series/rot_hist.mat')
poly = create_poly_struct(st_series);
colores=lines(size(rot_hist,2));

for i = 1:size(rot_hist,2)
    rs=[rot_hist{:,i}];

    subplot(4,2,2)
    plot(rs(1,:),rs(5,:),'color',colores(i,:))
    hold on
    title('Scale factor')

    subplot(4,2,4)
    plot(rs(1,:),rs(6,:),'color',colores(i,:))
    hold on
    title('X translation')

    subplot(4,2,6)
    plot(rs(1,:),rs(7,:),'color',colores(i,:))
    hold on
    title('Y translation')

    subplot(4,2,8)
    plot(rs(1,:),rs(8,:),'color',colores(i,:))
    hold on
    title('Z translation')
end

figure
clf
yticks = [];
for i = 1:size(stab_sites_i,1)
    t=[st_series(stab_sites_i(i)).epochs]';
    
    Ha=etm{stab_sites_i(i),2};
    if ~isempty(Ha)
        if ~isempty(min(Ha(Ha(:,2) ~= 0,1)))
            tp = t(t >= min(Ha(Ha(:,2) ~= 0,1)));

            post=plot(tp,repmat(i,[size(tp,1) 1]),'or','MarkerFaceColor','r','MarkerSize',3);

            hold on

            ti = t(t < min(Ha(Ha(:,2) ~= 0,1)));
            plot(ti,repmat(i,[size(ti,1) 1]),'ob','MarkerFaceColor','b','MarkerSize',3)
        else
            inter=plot(t,repmat(i,[size(t,1) 1]),'ob','MarkerFaceColor','b','MarkerSize',3);
        end
    else
        inter=plot(t,repmat(i,[size(t,1) 1]),'ob','MarkerFaceColor','b','MarkerSize',3);
    end
    
    yticks=[yticks; st_series(stab_sites_i(i)).stnm];
end

set(gca,'YTick',1:size(stab_sites_i,1),'YTickLabel',yticks)
set(gca,'FontSize',6,'FontName','courier')
grid on

legend([inter post],'Inter-seismic','Post-seismic')

st=poly.x(:,stab_sites_i);
numstab = size(st,2) - sum(isnan(st),2);

title(['Stabilization sites history: Max stab sites: ' num2str(max(numstab)) '; Min stab sites: ' num2str(min(numstab))],'FontSize',14)

[~,m]=max(numstab);
plot([poly.epochs(m) poly.epochs(m)],ylim,'g');

[~,m]=min(numstab);
plot([poly.epochs(m) poly.epochs(m)],ylim,'r');

