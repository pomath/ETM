clc
clear all

% load the stack solution
load('../ETM_files/st_series/stab_sites_i.mat');
load('../ETM_files/st_series/st_series.mat');
load('../ETM_files/st_series/etm.mat');

% load the station info
load('../ETM_files/st_series/st_info.mat');

% load the SINEX solution ref sites
[ref_sites,ref_sites_i, ref_sites_v] = load_ref_vel(st_series,etm,'../tables/SIR09P01.SNX');

% align stacked RF in velocity space
poly = create_poly_struct(st_series);

X = []; Y = []; Z = [];
vx = []; vy = []; vz = [];
t = nan(size(poly.epochs,1),size(poly.x,2));
% loop through the stations and get one coordinate for each one
for i = 1:size(st_series,2)
    if ismember(i,ref_sites_i)
        vx = [vx; etm{i,1}(1,2) - ref_sites_v(ismember(ref_sites_i,i),1)];
        vy = [vy; etm{i,1}(2,2) - ref_sites_v(ismember(ref_sites_i,i),2)];
        vz = [vz; etm{i,1}(3,2) - ref_sites_v(ismember(ref_sites_i,i),3)];
        index(i) = 1;
        % save the velocities to compare them later
        setm_vxyz(ismember(ref_sites_i,i),1) = etm{i,1}(1,2);
        setm_vxyz(ismember(ref_sites_i,i),2) = etm{i,1}(2,2);
        setm_vxyz(ismember(ref_sites_i,i),3) = etm{i,1}(3,2);
    else
        index(i) = 0;
    end
    % get a mean value for all stations, even if it doesn't have an etm
    X(1,i) = mean(st_series(i).x);
    Y(1,i) = mean(st_series(i).y);
    Z(1,i) = mean(st_series(i).z);
    t(ismember(poly.epochs,st_series(i).epochs),i) = st_series(i).epochs - min(st_series(i).epochs);
end

index = logical(index);

n = size(X,2);

% build the design matrix
Ax = [zeros(n,1) -Z' Y'];
Ay = [Z' zeros(n,1) -X'];
Az = [-Y' X' zeros(n,1)];

A = [Ax; Ay; Az];
Ai = [Ax; Ay; Az];

% remove the stations with no ETM
A([index'; index'; index'] == 0,:) = [];

rt = [poly.x poly.y poly.z];

% L has the velocities
L = [vx; vy; vz];
x = (A'*A)\(A'*L);

% remove the adjusted velocities
rt = rt - [t t t].*repmat((Ai*x)',size(t,1),1);

st_series_o = load_poly2series(st_series,rt,[]);

[etm, index, weights] = auto_fit_xyz(st_series_o,st_info,ref_sites_i);
st_series_o = update_st_series(st_series_o,index,weights);

for i = 1:size(st_series,2)
    if ismember(i,ref_sites_i)
        % to plot the final velocities
        etm_vxyz(ismember(ref_sites_i,i),1) = etm{i,1}(1,2);
        etm_vxyz(ismember(ref_sites_i,i),2) = etm{i,1}(2,2);
        etm_vxyz(ismember(ref_sites_i,i),3) = etm{i,1}(3,2);
    end
end
        
clf
hold on
lon=[st_series(ref_sites_i).lon]';
lon(lon > 180) = lon(lon > 180) - 360;
lat=[st_series(ref_sites_i).lat]';
stnm=char({st_series(ref_sites_i).stnm});
m_proj('miller','lat',[min(lat)-1 max(lat)+1],'long',[min(lon)-1 max(lon)+1])
m_coast();
m_plot(lon,lat,'or');
m_text(lon,lat,stnm,'FontSize',6);

% convert velocities to NEU
[ref_vneu(:,1),ref_vneu(:,2),ref_vneu(:,3)] = ct2lg(ref_sites_v(:,1),ref_sites_v(:,2),ref_sites_v(:,3),lat*pi/180,lon*pi/180);
[etm_vneu(:,1),etm_vneu(:,2),etm_vneu(:,3)] = ct2lg(etm_vxyz(:,1),etm_vxyz(:,2),etm_vxyz(:,3),lat*pi/180,lon*pi/180);
[setm_vneu(:,1),setm_vneu(:,2),setm_vneu(:,3)] = ct2lg(setm_vxyz(:,1),setm_vxyz(:,2),setm_vxyz(:,3),lat*pi/180,lon*pi/180);

m_quiver(lon,lat,ref_vneu(:,2),ref_vneu(:,1),'b')
m_quiver(lon,lat,setm_vneu(:,2),setm_vneu(:,1),'g')
m_quiver(lon,lat,etm_vneu(:,2),etm_vneu(:,1),'r')

axis equal
axis tight