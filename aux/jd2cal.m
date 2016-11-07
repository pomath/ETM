function [yr,mn,dy]=jd2cal(jd)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

a = fix(jd+0.5);
if a < 2299161
    c = a + 1524;
else
    b = fix( (a-1867216.25) / 36524.25 );
    c = a + b - fix(b/4) + 1525;
end
d = fix( (c-122.1)/365.25 );
e = fix(365.25*d);
f = fix( (c-e) / 30.6001 );
dy = c - e - fix(30.6001*f) + rem((jd+0.5),a);
mn = f - 1 - 12*fix(f/14);
yr = d - 4715 - fix( (7+mn)/10 );

end

