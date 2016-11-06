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
[poly,st_series] = remove_stn('LLNF', poly,st_series);
%[poly,st_series] = remove_stn('MCLA', poly,st_series);
%[poly,st_series] = remove_stn('CHLT', poly,st_series);

% load the station info
st_info = load_st_info('../tables/station.info');

load 'st_series'

X = poly.x;
Y = poly.y;
Z = poly.z;

Lr = [];
del_epochs = [];

for i = 2:size(poly.epochs,1)
    
    % build the design matrix with the coordinates of the current epoch
    % we want to find the rot and trans that produce the best fit between
    % the previous and current epoch
    
    Ax = [zeros(size(poly.x,2),1) Z(i,:)' -Y(i,:)' repmat([1 0 0],size(poly.x,2),1)];
    Ay = [-Z(i,:)' zeros(size(poly.x,2),1) X(i,:)' repmat([0 1 0],size(poly.x,2),1)];
    Az = [Y(i,:)' -X(i,:)' zeros(size(poly.x,2),1) repmat([0 0 1],size(poly.x,2),1)];

    % verify if there is a co-seismic jump between epoch i and i-1
    %Aj = diag(load_hsf_all(poly.epochs(i-1),poly.epochs(i),char(extractfield(st_series,'stnm'))));
    Aj = load_hsf_all(poly.epochs(i-1),poly.epochs(i),char(extractfield(st_series,'stnm')));
    % remove all columns with zeros
    % Aj(:,~any(Aj,1)) = [];
    A = [Ax; Ay; Az];
    
    % select the previous epoch
    L = [X(i-1,:)' - X(i,:)'; Y(i-1,:)' - Y(i,:)'; Z(i-1,:)' - Z(i,:)'];
    
    % if the result is not empty, attach to Ax Ay Az
    if sum(Aj) > 0
        % if it has a jump, remove it from the adjustment
        L(logical(Aj),:) = [];
        A(logical(Aj),:) = [];
        %A = [Ax Aj zeros(size(Aj,1),size(Aj,2)*3-size(Aj,2)); ...
        %     Ay zeros(size(Aj,1),size(Aj,2)) Aj zeros(size(Aj,1),size(Aj,2)); ...
        %     Az zeros(size(Aj,1),size(Aj,2)*3-size(Aj,2)) Aj];
    % else
        %A = [Ax; Ay; Az];
    end
    
    L(any(isnan(sum(A,2)),2),:)=[];
    A(any(isnan(sum(A,2)),2),:)=[];
    % now remove from A rows where L = NaN
    % this happens because of the 3*sigma filtered time serie
    % some data points get removed from the TS but not from the
    % polyhedra.
    A(isnan(L(:,1)),:)=[];
    L(isnan(L(:,1)),:)=[];
    
    % after everything is done, check that there are no zero cols left in A
    % this can happen when a jump is introduced but there is no data point
    % for the station (therefore, the row is removed)
    A(:,~any(A,1)) = [];
    
    if size(A,1) > size(A,2)
        x = (A'*A)\(A'*L);
        
        % make A again to include the NaN back
        A = [Ax; Ay; Az];
        % reload L to filter A again
        L = [X(i-1,:)' - X(i,:)'; Y(i-1,:)' - Y(i,:)'; Z(i-1,:)' - Z(i,:)'];
        % Put NaNs in the position of the outliers
        % A(isnan(L(:,1)),:)=NaN;

        % apply the rotation and translation
        Lr = [Lr; (A*x(1:6))' + [X(i,:) Y(i,:) Z(i,:)]];
    else
        % insufficient observation equations!
        % remove points from time series
        Lr = [Lr; nan(1,size(Lx,2).*3)];
        % save a vector of the epochs to delete
        del_epochs = [del_epochs; poly.epochs(i)];
    end
end

st_series_r = load_poly2series(st_series,Lr,del_epochs);

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
