% PLOT DMS PROFILES FROM GREEN EDGE
clear
clc

load gcms_greenedge_proc
load ALLSURF
rst = ALLSURF(:,strcmp('station',headALLSURF));

%% INITIAL SETTINGS

col.DMS = 46;
col.DMSPt = 49;
colors = brewermap(10,'Spectral');
bubsc.DMS = 40; % bubble plot scaling factor
bubsc.DMSPt = 4; % equal size is DMS:DMSPt = 0.1

%% PREPARE DATA

DATA = [data DMSOUT];
% dateaxstr = cellstr(datestr(dateax));

% Create station list
station = DATA(:,6);
diffstation = diff(station);
stnax = [station(diffstation~=0); station(end)];
stnaxstr = cellstr(num2str(stnax));
%
% % Preallocate station surface stats
% SURF = [stnax nan(length(stnax),1)];

% Define transects and longitudes
transects = {[400 403 406 409 413 412 418],...
    [507 506 512 515 519],...
    [600 603 605 615 612 604.5 608],...
    [703 700 707 713 716 719]};
% xmax.DMS = [10 30 30 10];
xmax.DMS = [24 75 46 12];
xmax.DMSPt = [260 600 305 200];
zmax = [62 72 80 80];
lat = [68 70 70.5 69.5];
lon = {[-62.42 -61.6 -60.82 -60 -58.93 -59.2 -57.77],...
    [-59.1 -58.66 -60.35 -61.25 -62.42],...
    [-64 -63.05 -62.42 -59.5 -60.42 -62.63 -61.62],...
    [-58.72 -57.87 -59.8 -61.58 -62.43 -63.23]};
fignames = {'400' '500' '600' '700'};
% ytext = 22*zmax/80;
ytext = -10*ones(1,length(transects));

%% LOOP FOR PLOTTING PROFILES AND TRANSECTS

for j = 1:length(lon) % loop on transects
    
    jlons = lon{j};
    jtrans = transects{j};
    jlat = lat(j);
    
    % Sort stations by longitude
    l_t = [jlons' jtrans']; l_t = sortrows(l_t,1);
    jlons = l_t(:,1); jtrans = l_t(:,2);
    
    % Define colors depending on transect length
    tl = length(jtrans);
    if ~mod(tl,2)
        jcolors = [colors(1:tl/2,:); colors((end-tl/2+1):end,:)];
    else
        jcolors = [colors(1:floor(tl/2),:); colors((end-floor(tl/2)):end,:)];
    end
    jcolors = flip(jcolors,1);
    
    % Initialize ALLSURF data match
    RTRANS = [];
    
    figure(j), clf
    set(gcf,'units','centimeters','position',[1 2 42 17.5]) % MBP screen
    
    for k = 1:length(jlons) % loop on transect stations
        
        kstn = jtrans(k);
        tmp = DATA(station==kstn,:);
        RTRANS = [RTRANS; nanmean(ALLSURF(rst==kstn,:),1)];
        
        if ~isempty(tmp)
            
            date = datenum(tmp(1,1:3)); date_s = datestr(date); % date
            d_m_y = regexp(date_s,'-','split');
            
            % Assign variables for plotting. Select samples according to
            % pre-screening applied (100 in tr 400 and 500, 5 µm Nitex afterwards)
            % Separate casts
            ps = tmp(:,15); % prescreen
            cc = tmp(:,7); % cast
            c1 = cc(1);
            ce = cc(end);
            if max(jtrans) < 600
                z = tmp(ps==100&cc==c1,13);
                x = tmp(ps==100&cc==c1,col.DMS); xsd = tmp(ps==100&cc==c1,col.DMS+1);
                x2 = tmp(ps==100&cc==c1,col.DMSPt);
                if ce ~= c1
                    zb = tmp(ps==100&cc==ce,13);
                    xb = tmp(ps==100&cc==ce,col.DMS); xsdb = tmp(ps==100&cc==ce,col.DMS+1);
                    xb2 = tmp(ps==100&cc==ce,col.DMSPt);
                end
            else
                z = tmp(ps==5&cc==c1,13);
                x = tmp(ps==5&cc==c1,col.DMS); xsd = tmp(ps~=100&cc==c1,col.DMS+1);
                x2 = tmp(ps==5&cc==c1,col.DMSPt);
                if ce ~= c1
                    zb = tmp(ps==5&cc==ce,13);
                    xb = tmp(ps==5&cc==ce,col.DMS); xsdb = tmp(ps~=100&cc==ce,col.DMS+1);
                    xb2 = tmp(ps==5&cc==ce,col.DMSPt);
                end
            end
            
            tmp1 = [z x xsd];
            tmp2 = [z x2];
            tmp1 = sortrows(tmp1(~isnan(x),:),1);
            tmp2 = sortrows(tmp2(~isnan(x2),:),1);
            z = tmp1(:,1); x = tmp1(:,2); xsd = tmp1(:,3);
            z2 = tmp2(:,1); x2 = tmp2(:,2);
            % Add second cast if present
            if ce ~= c1
                tmp1_b = [zb xb xsdb];
                tmp2_b = [zb xb2];
                tmp1_b = sortrows(tmp1_b(~isnan(xb),:),1);
                tmp2_b = sortrows(tmp2_b(~isnan(xb2),:),1);
                zb = tmp1_b(:,1); xb = tmp1_b(:,2); xsdb = tmp1_b(:,3);
                zb2 = tmp2_b(:,1); xb2 = tmp2_b(:,2);
            end
            
            % ///// Bubble plot /////
            subplot(2,3,4)
            % DMS
            hold on
            scatter(jlons(k)*ones(size(z)),-z,bubsc.DMS*(x),jcolors(k,:),'filled') % +0.1
            if ce ~= c1
                scatter(jlons(k)*ones(size(zb)),-zb,bubsc.DMS*(xb),jcolors(k,:),'filled') % +0.05
            end
            % DMSPt
            scatter(jlons(k)*ones(size(z2)),-z2,bubsc.DMSPt*(x2),'k','linewidth',2) % +0.1
            if ce ~= c1
                scatter(jlons(k)*ones(size(zb2)),-zb2,bubsc.DMSPt*(xb2),'k','linewidth',2) % +0.05
            end
            
            ax1 = gca;
            ax12 = axes('Position',get(ax1,'Position'),'xtick',jlons,'XTickLabel',jtrans,'XAxisLocation','top','XColor','k','YColor','k');
            linkaxes([ax1 ax12],'x');
            
            % ///// Profile plot, DMS /////
            subplot(2,3,5), hold on
            for ii = 1:length(z), plot([x(ii)-xsd(ii) x(ii)+xsd(ii)],[-z(ii) -z(ii)],'-','color',jcolors(k,:),'linewidth',2), end % errorbars
            plot(x,-z,'-o','color',jcolors(k,:),'markersize',8,'markerfacecolor',jcolors(k,:),'linewidth',2)
            if ce ~= c1
                for ii = 1:length(zb), plot([xb(ii)-xsdb(ii) xb(ii)+xsdb(ii)],[-zb(ii) -zb(ii)],'-','color',jcolors(k,:),'linewidth',2), end % errorbars
                plot(xb,-zb,'o','color',jcolors(k,:),'markersize',8,'markerfacecolor',jcolors(k,:),'linewidth',2)
            end
            
            % ///// Profile plot, DMSPt /////
            subplot(2,3,6), hold on
            plot(x2,-z2,'-o','color',jcolors(k,:),'markersize',8,...
                'markeredgecolor',jcolors(k,:),'markerfacecolor','w','linewidth',2)
            if ce ~= c1
                plot(xb2,-zb2,'o','color',jcolors(k,:),'markersize',8,...
                    'markeredgecolor',jcolors(k,:),'markerfacecolor','w','linewidth',2)
            end
            
        end
    end
    
    % Variables from Achim's dataset
    iso04 = RTRANS(:,strcmp('isolume_m_at_0.415_Einm2-d1',headALLSURF));
    mld = RTRANS(:,strcmp('mld03',headALLSURF)); % 'MixedLayerDepth_0.1gkg-1-criterion_m'
    hBD = RTRANS(:,strcmp('hBD_m',headALLSURF));
    no3c = RTRANS(:,strcmp('Nitracline_m',headALLSURF));
    ice = RTRANS(:,strcmp('IceConcentration_percent',headALLSURF));
    dbm = RTRANS(:,strcmp('dbm',headALLSURF));
    owd = RTRANS(:,strcmp('OWD',headALLSURF));
    
    subplot(2,3,1)
    % --------- If showing Ice with area, OWD in letters ---------
    %     area(jlons(~isnan(ice)), ice(~isnan(ice)),'facecolor',[.7 .7 .7])
    %     set(gca,'units','centimeters','position',[4 13 12 1.2],'fontsize',18,'fontname','arial',...
    %         'xtick',jlons,'xticklabel',[],'xaxislocation','bottom',...
    %         'ytick',[0 100])
    %     ylabel('Ice %')
    %     axis([min(jlons) max(jlons) 0 100])
    %     grid on, box off
    %     % Text
    %     for k = 1:length(jlons) % loop on transect stations
    %         % Text displaying St. number, and either
    %         %             % date
    %         %             text(jlons(k),ytext(j),sprintf('St. %s, %s-%s',num2str(jtrans(k)),d_m_y{1},d_m_y{2}),'color',jcolors(k,:),...
    %         %                 'fontweight','bold','fontsize',14,'fontname','arial','rotation',45)
    %         % or OWD
    %         text(jlons(k),ytext(j),sprintf('St. %s, OWD = %i',num2str(jtrans(k)),tmp_owd),'color',jcolors(k,:),...
    %             'fontweight','bold','fontsize',14,'fontname','arial','rotation',45)
    %     end
    
    % --------- If showing Ice with bars, OWD in points-lines ---------
    bar(jlons(~isnan(ice)),ice(~isnan(ice)),'facecolor',[.8 .8 .8],'edgecolor',[.8 .8 .8],...
        'barwidth',0.6,'linewidth',2); hold on
    lowd = [nanmin(owd) nanmax(owd)];
    towd = floor(lowd(1)/10)*10:10:ceil(lowd(2)/10)*10;
    sowd = 95*(owd-lowd(1))/diff(lowd); % scale to Ice axis
    plot(jlons(~isnan(sowd)),sowd(~isnan(sowd)),'o-k','markersize',10,...
        'markerfacecolor','w','linewidth',2)
    for k = 1:length(jlons)
        jk = jlons(k);
        ok = owd(k);
        sok = sowd(k);
        text(jk(~isnan(ok)),sok(~isnan(sok))+17,sprintf('%0.0f',ok),'fontsize',12)
        text(jlons(k),ytext(j),sprintf('%s',num2str(jtrans(k))),...
            'color',jcolors(k,:),'fontweight','bold','fontsize',14,'fontname','arial','rotation',-60)
    end
    text(max(jlons),max(sowd(~isnan(sowd)))-20,'Open Water Days','fontsize',16,'rotation',0)
    set(gca,'units','centimeters','position',[3 13.8 14 3],'fontsize',18,'fontname','arial',...
        'xtick',jlons,'xticklabel',[],'xaxislocation','bottom',...
        'ytick',[0 100])
    ylabel('Ice %')
    axis([min(jlons)-0.5 max(jlons)+0.5 0 100])
    grid off, box off
    
%     ax12 = gca; % current axes
%     ax13 = axes('Position',get(ax12,'Position'),...
%         'Xlim',get(ax12,'XLim'),'XAxisLocation','bottom','Color','none','XColor','k',...
%         'Ylim',lowd, 'YAxisLocation','right','YColor','k','ytick',towd);
%     plot(jlons(~isnan(owd)),owd(~isnan(owd)),'-k','linewidth',2,'Parent',ax13,'Color','k')
    
    subplot(2,3,4)
    set(gca,'units','centimeters','position',[3 2 14 10],'fontsize',18,'fontname','arial')
    axis([min(jlons)-0.5 max(jlons)+0.5 -zmax(j) 0])
    xlabel(ax1,'Longitude','fontweight','bold','fontsize',18,'fontname','arial')
    ylabel('Depth (m)','fontweight','bold','fontsize',18,'fontname','arial')
    box on
    ax1 = gca;
    % p5 = plot(jlons(~isnan(mld)), -mld(~isnan(mld)),':k','linewidth',2,'parent',ax1);
    p4 = plot(jlons(~isnan(hBD)), -hBD(~isnan(hBD)),'--k','linewidth',2,'parent',ax1);
    p3 = plot(jlons(~isnan(iso04)), -iso04(~isnan(iso04)),'--b','linewidth',2,'parent',ax1);
    p2 = plot(jlons(~isnan(no3c)), -no3c(~isnan(no3c)),'--m','linewidth',2,'parent',ax1);
    p1 = plot(jlons(~isnan(dbm)), -dbm(~isnan(dbm)),':k','linewidth',2,'parent',ax1);
    uistack(p1,'bottom'); uistack(p2,'bottom'); uistack(p3,'bottom'); uistack(p4,'bottom'); 
    % uistack(p5,'bottom')
    if j >= 1
        % l = legend('XLD','hBD','Isolume','NO3-cline','DBM'); % MLD_{0.1}
        l = legend('hBD','Isolume','NO_{3}-cline','DBM'); % MLD_{0.1}
        set(l,'fontsize',16)
    end
    
    subplot(2,3,5)
    set(gca,'units','centimeters','position',[19 2 10 10],'fontsize',18,'fontname','arial',...
        'yticklabel','','xaxislocation','top')
    axis([0 xmax.DMS(j) -zmax(j) 0])
    xlabel(sprintf('%s\n','DMS (nM)'),'fontweight','bold','fontsize',18,'fontname','arial')
    box on
    
    subplot(2,3,6)
    set(gca,'units','centimeters','position',[30 2 10 10],'fontsize',18,'fontname','arial',...
        'yticklabel','','xaxislocation','top')
    axis([0 xmax.DMSPt(j) -zmax(j) 0])
    xlabel(sprintf('%s\n','DMSPt (nM)'),'fontweight','bold','fontsize',18,'fontname','arial')
    box on
    
end

