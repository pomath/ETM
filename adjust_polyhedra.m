clc
clear all
dy_series = load_daily();

% load the time series
[poly,st_series] = load_polyhedra();

% remove stations with few obs
rm_stn = [];
for i = 1:size(st_series,2)
    if size(st_series(i).epochs,2) < 50 | (max(st_series(2).epochs) - min(st_series(2).epochs)) < 1
        rm_stn = [rm_stn; st_series(i).stnm];
    end
end

for i = 1:size(rm_stn,1)
    [poly,st_series] = remove_stn(rm_stn(i,:), poly,st_series);
end

[poly,st_series] = remove_stn('CORD', poly,st_series);
[poly,st_series] = remove_stn('ACPM', poly,st_series);
[poly,st_series] = remove_stn('CFAG', poly,st_series);
[poly,st_series] = remove_stn('EMAK', poly,st_series);
[poly,st_series] = remove_stn('EPRF', poly,st_series);
[poly,st_series] = remove_stn('GAS1', poly,st_series);
[poly,st_series] = remove_stn('IGM0', poly,st_series);
[poly,st_series] = remove_stn('UCOR', poly,st_series);
[poly,st_series] = remove_stn('BUE2', poly,st_series);
[poly,st_series] = remove_stn('FEDE', poly,st_series);
[poly,st_series] = remove_stn('ARCO', poly,st_series);
[poly,st_series] = remove_stn('LNDS', poly,st_series);
[poly,st_series] = remove_stn('ELA2', poly,st_series);
[poly,st_series] = remove_stn('RMLS', poly,st_series);
[poly,st_series] = remove_stn('JVGO', poly,st_series);
[poly,st_series] = remove_stn('DORE', poly,st_series);
[poly,st_series] = remove_stn('MZAL', poly,st_series);
[poly,st_series] = remove_stn('MZGA', poly,st_series);
[poly,st_series] = remove_stn('PATA', poly,st_series);
[poly,st_series] = remove_stn('PDE3', poly,st_series);
[poly,st_series] = remove_stn('PECL', poly,st_series);
[poly,st_series] = remove_stn('PHDP', poly,st_series);
[poly,st_series] = remove_stn('RECO', poly,st_series);
[poly,st_series] = remove_stn('RSCL', poly,st_series);
[poly,st_series] = remove_stn('RUFI', poly,st_series);
[poly,st_series] = remove_stn('RWSN', poly,st_series);
[poly,st_series] = remove_stn('SCLA', poly,st_series);
[poly,st_series] = remove_stn('SOLD', poly,st_series);
[poly,st_series] = remove_stn('SURY', poly,st_series);
[poly,st_series] = remove_stn('TAVA', poly,st_series);
[poly,st_series] = remove_stn('TILC', poly,st_series);
[poly,st_series] = remove_stn('TMCO', poly,st_series);
[poly,st_series] = remove_stn('VIMA', poly,st_series);
[poly,st_series] = remove_stn('i5MA', poly,st_series);
[poly,st_series] = remove_stn('iARO', poly,st_series);
[poly,st_series] = remove_stn('EISL', poly,st_series);
[poly,st_series] = remove_stn('ESCA', poly,st_series);
[poly,st_series] = remove_stn('CMPN', poly,st_series);
[poly,st_series] = remove_stn('LLFN', poly,st_series);
[poly,st_series] = remove_stn('YEMA', poly,st_series);
[poly,st_series] = remove_stn('DUAO', poly,st_series);
[poly,st_series] = remove_stn('CAUQ', poly,st_series);
[poly,st_series] = remove_stn('NIHU', poly,st_series);

% load the station info
st_info = load_st_info('../tables/station.info');

% fit etms and filter large outliers
[~,~,~,st_series_u] = auto_fit_xyz(st_series,st_info,poly.epochs,0,1);

st_series = st_series_u;

% refit the etms with no outliers
[Lx,Ly,Lz,~] = auto_fit_xyz(st_series,st_info,poly.epochs,0,0);

X = poly.x;
Y = poly.y;
Z = poly.z;

ref_stn = ['BRAZ';'OHI2'; 'AUTF'; 'PALM'; 'NAUS'; 'PALM'; 'FALK'; 'PARC'; 'RIO2'; 'GLPS'];
ref_index = find(ismember(extractfield(st_series,'stnm'),ref_stn));

for j = 1:5
    % matrix to put the results
    Lr = [];
    del_epochs = [];
    % build the design matrix
    for i = 1:size(poly.epochs,1)
        % loop through the epochs

        Ax = [zeros(size(poly.x,2),1) Z(i,:)' -Y(i,:)' repmat([1 0 0],size(poly.x,2),1)];
        Ay = [-Z(i,:)' zeros(size(poly.x,2),1) X(i,:)' repmat([0 1 0],size(poly.x,2),1)];
        Az = [Y(i,:)' -X(i,:)' zeros(size(poly.x,2),1) repmat([0 0 1],size(poly.x,2),1)];

        L = [Lx(i,:)' - X(i,:)'; Ly(i,:)' - Y(i,:)'; Lz(i,:)' - Z(i,:)'];

        A = [Ax; Ay; Az];
        
        p = [poly.px(i,:)'; poly.py(i,:)'; poly.pz(i,:)'];
        
        % remove rows with zeros (those that don't have observations)
        % L(all(A(:,1:9)==0,2),:)=[];
        p(any(isnan(sum(A,2)),2),:)=[];
        L(any(isnan(sum(A,2)),2),:)=[];
        % A(all(A(:,1:9)==0,2),:)=[];
        A(any(isnan(sum(A,2)),2),:)=[];
        % now remove from A rows where L = NaN
        % this happens because of the 3*sigma filtered time serie
        % some data points get removed from the TS but not from the
        % polyhedra.
        A(isnan(L(:,1)),:)=[];
        p(isnan(L(:,1)),:)=[];
        L(isnan(L(:,1)),:)=[];
        
        if size(A,1) > size(A,2)
            
            P = diag(1./(p.^2));
            x = (A'*A)\(A'*L);

            % make A again to include the NaN back
            A = [Ax; Ay; Az];
            % reload L to filter A again
            L = [Lx(i,:)'-X(i,:)'; Ly(i,:)'-Y(i,:)'; Lz(i,:)'-Z(i,:)'];
            % Put NaNs in the position of the outliers
            A(isnan(L(:,1)),:)=NaN;

            % apply the rotation and translation
            Lr = [Lr; (A*x)' + [X(i,:) Y(i,:) Z(i,:)]];
        
            % the new X Y Z coordinates to perform the next adjustment are the rotated ones
            X(i,:) = (Ax*x)' + X(i,:);
            Y(i,:) = (Ay*x)' + Y(i,:);
            Z(i,:) = (Az*x)' + Z(i,:);
        else
            % insufficient observation equations!
            % remove points from time series
            Lr = [Lr; nan(1,size(Lx,2).*3)];
            % save a vector of the epochs to delete
            del_epochs = [del_epochs; poly.epochs(i)];
        end
    end

    % rebuild the st_series from the Lr matrix
    st_series_r = load_poly2series(st_series,Lr,del_epochs);
    
    % refit without common mode
    [Lx,Ly,Lz,~,osc_param] = auto_fit_xyz(st_series_r,st_info,poly.epochs,0,0);
    
    disp(['Run ' num2str(j) ' of 5'])
end

% save the ts
st_series_bak = st_series_r;

% remove common mode
st_series_r = remove_cmm(st_series_r,poly,osc_param);
    
% convert the results in st_series_r to NEU
for i = 1:size(st_series_r,2)
    dx = st_series_r(i).x' - mean(st_series_r(i).x);
    dy = st_series_r(i).y' - mean(st_series_r(i).y);
    dz = st_series_r(i).z' - mean(st_series_r(i).z);
    
    [n,e,u]=ct2lg(dx,dy,dz,st_series_r(i).lat*pi/180,st_series_r(i).lon*pi/180);
    
    st_series_r(i).n = n';
    st_series_r(i).e = e';
    st_series_r(i).d = u';
end
