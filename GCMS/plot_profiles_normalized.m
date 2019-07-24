% PLOT DMS PROFILES FROM GREEN EDGE
clear
clc
close all

load gcms_greenedge_proc
load ~/Desktop/GreenEdge/Irradiance/CASTS_ICE_matched_25-Apr-2017.mat

%% INITIAL SETTINGS

column = 45;
ncolors = 20;
% colors = brewermap(ncolors,'dark2');
colors = brewermap(ncolors,'spectral');
% colors = flip(brewermap(ncolors,'spectral'));

%% PREPARE DATA

DATA = [data DMSOUT];
% dateaxstr = cellstr(datestr(dateax));

% Create station list
station = DATA(:,6);
diffstation = diff(station);
stnax = [station(diffstation~=0); station(end)];
stnaxstr = cellstr(num2str(stnax));

% Define transects and longitudes
transects = {[400 403 409 413 412 418],...
    [507 506 512 515 519],...
    [600 603 605 615 612 604.5 608],...
    [703 700 707 713 716 719]};
% xmax2 = [10 30 30 10];
xmax2 = [24 75 46 12];
zmax = [62 72 100 80];
lat = [68 70 70.5 69.5];
lon = {[-62.42 -61.6 -60 -58.93 -59.2 -57.77],...
    [-59.1 -58.66 -60.35 -61.25 -62.42],...
    [-64 -63.05 -62.42 -59.5 -60.42 -62.63 -61.62],...
    [-58.72 -57.87 -59.8 -61.58 -62.43 -63.23]};
fignames = {'400' '500' '600' '700'};


%% LOOP FOR PLOTTING PROFILES AND TRANSECTS

figure(245), clf, hold on
set(gcf,'units','centimeters','position',[5 10 15 15])
count = 0;

for j = 1:length(lon) % loop on transects
    
    jlons = lon{j};
    jtrans = transects{j};
    jlat = lat(j);
    
    % Sort stations by longitude
    l_t = [jlons' jtrans']; l_t = sortrows(l_t,1);
    jlons = l_t(:,1); jtrans = l_t(:,2);
    
    for k = 1:length(jlons) % loop on transect stations
        
        kstn = jtrans(k);
        tmp = DATA(station==kstn,:);
        tmp(isnan(tmp(:,13)) & isnan(tmp(:,45)),:) = [];
        
        % Match cast to ICE
        tmp2 = OUT(OUT(:,8)==mode(tmp(:,7)),12:14);
        imean = nanmean(nanmean(tmp2))*ncolors;
        imean(imean<1) = 1;
        ccolor = interp1((1:ncolors),colors,imean);
        
        if size(tmp,1) > 2 && nanmax(tmp(:,13))>=35
            count = count + 1;
            z = tmp(:,13);
            x = tmp(:,45);
            % x = x/nanmean(x);
            % x = x/nanmedian(x);
            x = x/geomean(x);
            icolor = mod(count,ncolors);
            if mod(count,ncolors)==0, icolor = ncolors; end
            plot(x,-z,'.-','color',ccolor,'linewidth',1,'markersize',20), hold on, scatter(x,-z,40,ccolor,'filled')
            % plot(x,-z,'.-','color',ccolor,'linewidth',1,'markersize',20)
            % plot(x,-z,'.-','color',colors(icolor,:),'markersize',12,'linewidth',2)
            % plot(x,-z,'.-','color',[0 .5 1],'markersize',8)
        end
    end
end

box on
set(gca,'xaxislocation','top','xscale','lin','ytick',-80:20:0,'tickdir','out')
title('Geometric-mean-normalized profiles','fontsize',18)
xlabel('','fontsize',18)
ylim([-80 0])
ylabel('Depth (m)','fontsize',18)
set(gca,'fontsize',14)
colormap(colors)
hc = colorbar;
set(hc,'location','southoutside','xticklabel',{'.2' '.4' '.6' '.8' '1'},'fontsize',14)
xlabel(hc,'Ice conc.','fontsize',18)

