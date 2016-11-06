function [poly,st_series] = remove_stn(stnm, poly,st_series)
    index = structfind(st_series,'stnm',stnm);
    st_series(:,index) = [];
    poly.x(:,index) = [];
    poly.y(:,index) = [];
    poly.z(:,index) = [];
    poly.px(:,index) = [];
    poly.py(:,index) = [];
    poly.pz(:,index) = [];
end