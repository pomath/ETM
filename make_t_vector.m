function ts = make_t_vector(mint, maxt)
    [y1,m1,d1]=jd2cal(yr2jd(mint));
    [y2,m2,d2]=jd2cal(yr2jd(maxt));
    ts = datenum(y1, m1, floor(d1)):datenum(y2, m2, floor(d2));
    ts = decyear(ts');
end