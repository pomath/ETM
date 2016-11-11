function plot_station(stnm, i, index, lat, lon, t, Ha, Cx, Cy, Cz, Lx, Ly, Lz)

    display(['Plotting ' stnm '...'])

    ts = (min(t):1/365:max(t))';
    % build the A matrix for a complete time series
    [A,~] = load_hsf3(stnm, ts, false, Ha);

    % get the A matrix using t instead of ts to calculate the STD
    [Ai,~] = load_hsf3(stnm, t(index), false, Ha);

    sigmaN = std(Lx(index) - Ai*Cx)*1000;
    sigmaE = std(Ly(index) - Ai*Cy)*1000;
    sigmaU = std(Lz(index) - Ai*Cz)*1000;
    
    % remove rows of zeros
    r = ~all(A == 0,2);
    A = A(r,:);
    ts = ts(r);
    
    h = clf;
    subplot(3,1,1)

    for comp = 1:3
        subplot(3,1,comp)
        switch comp
            case 1
                L = Lx; C = Cx;
                label = 'N';
            case 2
                L = Ly; C = Cy;
                label = 'E';
            case 3
                L = Lz; C = Cz;
                label = 'U';
        end

        hold on

        plot(t,L-C(1),'ob','MarkerSize', 2,'MarkerFaceColor','b')
        % plot outliers in cyan
        plot(t(~index),L(~index)-C(1),'oc','MarkerSize', 2,'MarkerFaceColor','c')

        M = A*C;
        plot(ts,M-C(1),'r','LineWidth',0.5);

        axis tight
        set(gca,'FontSize',10);

        % plot the Ha
        for j = 1:size(Ha,1)
            if Ha(j,1) > min(t) && Ha(j,1) < max(t)
                if Ha(j,2) == 0
                    color_l = ':b';
                else
                    color_l = ':r';
                end

                plot([Ha(j,1) Ha(j,1)], ylim,color_l,'LineWidth',1)
            end
        end
        grid on

        ylabel([label ' [m]'])

    end

    subplot(3,1,1)

    title({[stnm '(' num2str(i) '): \sigma_N = ' num2str(sigmaN) ' mm; \sigma_E = ' num2str(sigmaE) ' mm; \sigma_U = ' num2str(sigmaU) ' mm']; ['Outliers detected: ' num2str(sum(index == 0)) ' of ' num2str(size(index,1)) ' data points']},'FontSize',10)        

    % k = waitforbuttonpress;
    print(h,'-dpng',['etm_fless/' stnm '_NEU.png'], '-r300')
end
