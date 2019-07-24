% COMPARE DIFFERENT INTERPOLATION APPROACHES TO OBTAIN PAR AT NOON
clc, clear

load ~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2016-ICECAMP_cyber.mat

tvec = datevec(par_ts.mtimeLOC);
doy = floor(yearday(tvec(:,3),tvec(:,2),tvec(:,1),tvec(:,4),tvec(:,5),tvec(:,6)));
doylist = doy(tvec(:,4)==0);
ndays = length(doylist);

OUT = nan(ndays-1,6);

lon = -63.78953;
official_diff = -4; % UTC-4
real_diff = lon/15;
local_noon = 12 + (official_diff - real_diff);
local_noon = (round(local_noon*10))/10;

%% For each day...
for id = 1:ndays-1
    
    % Daily time series at 1h resolution
    x1 = tvec(doy == doylist(id),4);
    y1 = par_ts.data(doy == doylist(id));
    x1 = [x1; x1(end)+1];
    y1 = [y1; y1(end)];
    
    % Daily time series at 3h resolution (simulating satellite-RT output)
    i3 = ones(3,1)*(1:8);
    i3 = i3(:);
    x3 = (0:3:21)';
    y3 = nan(size(x3));
    for ip = 1:8
        y3(ip) = nanmean(y1(i3 == ip));
    end
    x3 = [x3; x3(end)+3];
    y3 = [y3; y3(end)];
    
    
    % Interpolate 1h time series
    xi = 9:0.1:15; % x for interpolation
    y1c = interp1(x1+0.5,y1,xi,'pchip');
    y1s = interp1(x1+0.5,y1,xi,'spline');
    y3c = interp1(x3+1.5,y3,xi,'pchip');
    y3s = interp1(x3+1.5,y3,xi,'spline');
    
    % Fill output matrix
    inoon = local_noon == xi;
    OUT(id,1) = y1(x1 == floor(local_noon));
    OUT(id,2) = y3(x3/3 == floor(local_noon/3));
    OUT(id,3:6) = [y1c(inoon) y1s(inoon) y3c(inoon) y3s(inoon)];
    
    if ~mod(id,7)
        h = figure(doylist(id));
        set(gcf,'units','centimeters','position',[5 5 40 40])
        stairs(x1,y1,'color','b','linewidth',2), hold on
        stairs(x3,y3,'color','r','linewidth',2)
        plot(xi,y1c,'-b','linewidth',1)
        plot(xi,y1s,'--b','linewidth',1)
        plot(xi,y3c,'-r','linewidth',1)
        plot(xi,y3s,'--r','linewidth',1)
        legend('1h data','3h data','1h PCHIP','1h spline','3h PCHIP','3h spline','fontsize',10)
        % Markers for noontime value
        plot(local_noon,OUT(id,1),'xb','markersize',8,'linewidth',2)
        plot(local_noon,OUT(id,2),'xr','markersize',8,'linewidth',2)
        plot(local_noon,OUT(id,3),'ob','markersize',8,'linewidth',2)
        plot(local_noon,OUT(id,4),'sb','markersize',8,'linewidth',2)
        plot(local_noon,OUT(id,5),'or','markersize',8,'linewidth',2)
        plot(local_noon,OUT(id,6),'sr','markersize',8,'linewidth',2)
        set(gca,'xlim',[0 24],'xtick',0:3:24,'fontsize',12)
        xlabel('Local hour'), ylabel('PAR, µmol photons m^{-2} s^{-1}')
        print(sprintf('interpPAR_DOY%03i.png',doylist(id)),'-dpng')
        close(h)
    end
    
end

%% Scatterplots comparing different methods

x = OUT(:,1);
X = [ones(size(x)) x];
methods = {'3h_data' '1h_PCHIP' '1h_spline' '3h_PCHIP' '3h_spline'};

for ij = 1:5
    
    y = OUT(:,ij+1);
    M = X\y;
    xp = sortrows(x);
    XP = [ones(size(xp)) xp];
    yp = XP*M;
    
    Xbis = [ones(size(y)) y];
    Mbis = Xbis\x;
    
    h = figure(ij);
    set(gcf,'units','centimeters','position',[5 5 20 20])
    plot([0 nanmax(xp)],[0 nanmax(xp)],'-k'), hold on
    plot(x,y,'ok')
    plot(xp,yp,'-m')
    axis([0 nanmax(xp) 0 nanmax(xp)])
    title(sprintf('y = %0.03fx + %0.02f\nx = %0.03fy + %0.02f',M(2),M(1),Mbis(2),Mbis(1)))
    xlabel('1h PAR')
    ylabel(methods{ij},'interpreter','none')
    print(sprintf('compare_interpPAR_%s.png',methods{ij}),'-dpng')
    close(h)
    
end

