function [doy,fyear]= date2doy(y,m,d,h)
%date2doy date to day of year (DOY)
%
% USAGE:       doy= date2doy(y,m,d)
%
%             [doy,fyear]= date2doy(y,m,d,h)
%
% INPUT
%      y       year (integer)
%      m       month (integer)
%      d       day of month (integer)
%      h       hour of day (real number in range 0 - 24)
%
% OUTPUT
%    doy       day of year (integer)
%  fyear       year and fractional part of year
%              
% Note: fyear is fractional year under the assumption that
% the year has 365 days, except for leap years with 366 days,
% and therefore constitutes a non-uniform measure of duration.
%
% See also function doy2date

lday=  [31 59 90 120 151 181 212 243 273 304 334 365]; 
if rem(y,4)==0     % it is a leap year
    lday(2:12)=lday(2:12)+1;
end

if m==1
    doy=d;
else
    doy=lday(m-1)+d;
end
if nargin<4
    h=0;
end
    
if nargout==2
    fyear=yyddd2fyear(y,doy,h);
end

        