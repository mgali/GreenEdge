% PLOT DMS PROFILES FROM GREEN EDGE
clear
clc

load gcms_greenedge_proc

path_rand = '~/Desktop/GreenEdge/Achim_Randelhoff/data/Randelhoff-et-al_Ice-edge-paper_per-station_v0.1.csv';
RAND = csvread(path_rand, 1, 0);
headRAND = {'station' 'isolume_m_at_0.1_Einm-2d-1' 'isolume_m_at_0.415_Einm2-d1'...
    'IceConcentration_percent' 'MixedLayerDepth_0.1gkg-1-criterion_m' 'hBD_m'...
    'Nitracline_m' 'AW-ArW_clustering_coefficient' 'OWD' 'PAR_at_3m_Einm-2d-1'};
rst = RAND(:,strcmp('station',headRAND));

%% INITIAL SETTINGS

col.DMSPt = 49;
colors = brewermap(10,'Spectral');
zsurf = 5.5;
bubsc.DMSPt = 5; % bubble plot scaling factor

%% PREPARE DATA

DATA = [data DMSOUT ones(size(DMSOUT,1),1)];
% dateaxstr = cellstr(datestr(dateax));

% Create station list
station = DATA(:,6);
diffstation = diff(station);
stnax = [station(diffstation~=0); station(end)];
stnaxstr = cellstr(num2str(stnax));

% Preallocate station surface stats
SURF = [stnax nan(length(stnax),1)];

% Define transects and longitudes
transects = {[400 403 406 409 413 412 418],...
    [507 506 512 515 519],...
    [600 603 605 615 612 604.5 608],...
    [703 700 707 713 716 719]};
% xmax.DMSPt = [10 30 30 10];
xmax.DMSPt = [260 600 200 200];
zmax = [62 72 100 80];
lat = [68 70 70.5 69.5];
lon = {[-62.42 -61.6 -60.82 -60 -58.93 -59.2 -57.77],...
    [-59.1 -58.66 -60.35 -61.25 -62.42],...
    [-64 -63.05 -62.42 -59.5 -60.42 -62.63 -61.62],...
    [-58.72 -57.87 -59.8 -61.58 -62.43 -63.23]};
fignames = {'400' '500' '600' '700'};


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
    
    % Initialize Achim Randelhoff's data match
    RTRANS = [];
    
    figure(j), clf
    % set(gcf,'units','centimeters','position',[5 10 50 27]) % large screen
    set(gcf,'units','centimeters','position',[1 2 30 16]) % MBP screen
    
    for k = 1:length(jlons) % loop on transect stations
        
        kstn = jtrans(k);
        RTRANS = [RTRANS; RAND(rst==kstn,:)];
        tmp = DATA(station==kstn,:);
        
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
                z = tmp(ps==100&cc==c1,13); x = tmp(ps==100&cc==c1,col.DMSPt); xsd = tmp(ps==100&cc==c1,col.DMSPt+1); scm = tmp(ps==100&cc==c1,30);
                if ce ~= c1
                    zb = tmp(ps==100&cc==ce,13); xb = tmp(ps==100&cc==ce,col.DMSPt); xsdb = tmp(ps==100&cc==ce,col.DMSPt+1);  scmb = tmp(ps==100&cc==ce,30);
                end
            else
                z = tmp(ps==5&cc==c1,13); x = tmp(ps==5&cc==c1,col.DMSPt); xsd = tmp(ps~=100&cc==c1,col.DMSPt+1);  scm = tmp(ps==5&cc==c1,30);
                if ce ~= c1
                    zb = tmp(ps==5&cc==ce,13); xb = tmp(ps==5&cc==ce,col.DMSPt); xsdb = tmp(ps~=100&cc==ce,col.DMSPt+1); scmb = tmp(ps==5&cc==ce,30);
                end
            end
            
            tmp2 = [z x xsd scm];
            tmp2(isnan(x),:) = [];
            tmp2 = sortrows(tmp2,1);
            z = tmp2(:,1); x = tmp2(:,2); xsd = tmp2(:,3); scm = tmp2(:,4);
            % Add second cast if present
            if ce ~= c1
                tmp2_b = [zb xb xsdb scmb];
                tmp2_b(isnan(xb),:) = [];
                tmp2_b = sortrows(tmp2_b,1);
                zb = tmp2_b(:,1); xb = tmp2_b(:,2); xsdb = tmp2_b(:,3);  scmb = tmp2_b(:,4);
            end
            
            % ///// Bubble plot /////
            subplot(1,2,1)
            scatter(jlons(k)*ones(size(z)),-z,bubsc.DMSPt*(x+0.1),jcolors(k,:),'filled'), hold on
            if sum(scm)>0, scatter(jlons(k),-z(scm==1),bubsc.DMSPt*(x(scm==1)+0.5),'k','linewidth',3), end
            if ce ~= c1
                scatter(jlons(k)*ones(size(zb)),-zb,bubsc.DMSPt*(xb+0.05),jcolors(k,:),'filled')
                if sum(scmb)>0, scatter(jlons(k),-zb(scmb==1),bubsc.DMSPt*(xb(scmb==1)+0.5),'k','linewidth',3), end
            end
            text(jlons(k),4,sprintf('%s, %s-%s',num2str(jtrans(k)),d_m_y{1},d_m_y{2}),'color',jcolors(k,:),...
                'fontweight','bold','fontsize',18,'fontname','arial','rotation',45)
            
            ax1 = gca;
            ax2 = axes('Position',get(ax1,'Position'),'xtick',jlons,'XTickLabel',jtrans,'XAxisLocation','top','XColor','k','YColor','k');
            linkaxes([ax1 ax2],'x');
            
            % ///// Profile plot /////
            subplot(1,2,2), hold on
            for ii = 1:length(z), plot([x(ii)-xsd(ii) x(ii)+xsd(ii)],[-z(ii) -z(ii)],'-','color',jcolors(k,:),'linewidth',2), end % errorbars
            plot(x,-z,'-o','color',jcolors(k,:),'markersize',8,'markerfacecolor',jcolors(k,:),'linewidth',2)
            if sum(scm)>0, plot(x(scm==1),-z(scm==1),'ok','markersize',11,'linewidth',3), end
            if ce ~= c1
                for ii = 1:length(zb), plot([xb(ii)-xsdb(ii) xb(ii)+xsdb(ii)],[-zb(ii) -zb(ii)],'-','color',jcolors(k,:),'linewidth',2), end % errorbars
                plot(xb,-zb,'o','color',jcolors(k,:),'markersize',8,'markerfacecolor',jcolors(k,:),'linewidth',2)
                if sum(scmb)>0, plot(xb(scmb==1),-zb(scmb==1),'ok','markersize',11,'linewidth',3), end
            end
            
%             % Station surface data
%             if ce == c1, SURF(stnax==jtrans(k),2) = nanmean(x(z<=zsurf));
%             elseif ce ~= c1, SURF(stnax==jtrans(k),2) = nanmean([x(z<=zsurf); xb(zb<=zsurf)]);
%             end
            
        end
    end
    
    subplot(1,2,1)
    % set(gca,'units','centimeters','position',[3 2 30 20],'fontsize',20,'fontname','arial') % large screen
    set(gca,'units','centimeters','position',[3 2 14 10],'fontsize',18,'fontname','arial') % MBP screen %
    axis([min(jlons)-0.5 max(jlons)+0.5 -zmax(j) 0])
    xlabel(ax1,'Longitude','fontweight','bold','fontsize',18,'fontname','arial')
    ylabel('Depth (m)','fontweight','bold','fontsize',18,'fontname','arial')
    box on
    
    iso04 = RTRANS(:,strcmp('isolume_m_at_0.415_Einm2-d1',headRAND));
    mld01 = RTRANS(:,strcmp('MixedLayerDepth_0.1gkg-1-criterion_m',headRAND));
    hBD = RTRANS(:,strcmp('hBD_m',headRAND));
    no3c = RTRANS(:,strcmp('Nitracline_m',headRAND));
    p1 = plot(jlons(~isnan(mld01)), -mld01(~isnan(mld01)),':k','linewidth',2,'parent',ax1);
    p2 = plot(jlons(~isnan(hBD)), -hBD(~isnan(hBD)),'--k','linewidth',2,'parent',ax1);
    p3 = plot(jlons(~isnan(iso04)), -iso04(~isnan(iso04)),'--b','linewidth',2,'parent',ax1);
    p4 = plot(jlons(~isnan(no3c)), -no3c(~isnan(no3c)),'--m','linewidth',2,'parent',ax1);
    uistack(p1,'bottom'); uistack(p2,'bottom'); uistack(p3,'bottom'); uistack(p4,'bottom')
    legend('NO3cli','Iso_{0.415}','hBD','MLD_{0.1}')
    
    subplot(1,2,2)
    % set(gca,'units','centimeters','position',[30 2 9 20],'fontsize',20,'fontname','arial',...
    set(gca,'units','centimeters','position',[19 2 10 10],'fontsize',18,'fontname','arial',...
        'yticklabel','','xaxislocation','top') % 'position',[22 2 7 12],
    axis([0 xmax.DMSPt(j) -zmax(j) 0])
    xlabel('DMSPt (nM)','fontweight','bold','fontsize',18,'fontname','arial')
    box on
    
    % % Automatic saving: cuts figures, regardless of format
    % print(j,sprintf('~/Desktop/GreenEdge/GCMS/plots_profile_v201901/transect_%s_DMS.png',fignames{j}),'-dpng');
    
end

%% Calculate stats separately for under-ice, ice margin stations

% uice = ([403 407 409 413 512 515 519 600 603 604.5 605 719])';
% uice = [uice nan(size(uice))];
% margin = ([415 418 506 507 608 612 615 700 703])';
% margin = [margin nan(size(margin))];
% 
% for j = 1:length(uice)
%     ff = SURF(:,1)==uice(j);
%     if sum(ff)
%         uice(j,2) = SURF(ff,2);
%     end
% end
% for j = 1:length(margin)
%     ff = SURF(:,1)==margin(j);
%     if sum(ff)
%         margin(j,2) = SURF(ff,2);
%     end
% end

