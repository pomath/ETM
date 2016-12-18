function plot_base_line(st_series,poly,stnm1,stnm2)
    clf
    site1=find(strcmp({st_series.stnm}, stnm1)==1);
    site2=find(strcmp({st_series.stnm}, stnm2)==1);

    refx = nanmean(poly.x(:,site1)-poly.x(:,site2));
    refy = nanmean(poly.y(:,site1)-poly.y(:,site2));
    refz = nanmean(poly.z(:,site1)-poly.z(:,site2));

    dx = poly.x(:,site1)-poly.x(:,site2)-refx;
    dy = poly.y(:,site1)-poly.y(:,site2)-refy;
    dz = poly.z(:,site1)-poly.z(:,site2)-refz;

    subplot(4,1,1)
    plot(poly.epochs,dx,'o-b')
    ylim([-nanstd(dx)*3 nanstd(dx)*3])
    grid on
    title(['Componenete X (' stnm1 ' - ' stnm2 ')'])
    
    subplot(4,1,2)
    plot(poly.epochs,dy,'o-g')
    ylim([-nanstd(dy)*3 nanstd(dy)*3])
    grid on
    title(['Componenete Y (' stnm1 ' - ' stnm2 ')'])
    
    subplot(4,1,3)
    plot(poly.epochs,dz,'o-r')
    ylim([-nanstd(dz)*3 nanstd(dz)*3])
    grid on
    title(['Componenete Z (' stnm1 ' - ' stnm2 ')'])
    
    subplot(4,1,4)
    D=sqrt(dx.^2 + dy.^2 + dz.^2);
    plot(poly.epochs,D-nanmean(D),'o-k')
    ylim([-nanstd(D)*3 nanstd(D)*3])
    grid on
    title(['Vector (' stnm1 ' - ' stnm2 ')'])
