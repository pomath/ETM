function yr=jd2yr(jd);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[iyr,mn,yr] = jd2cal(jd);
jd0 = cal2jd(iyr,1,1);
jd1 = cal2jd(iyr+1,1,1);
yr = iyr + (jd-jd0)./(jd1-jd0);

end

