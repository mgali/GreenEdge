% PLOT FDMS AND VERTICAL INTEGRALS, GREEN EDGE
clear
clc

% % Get here coordinates and doys for interpolation if needed
% load samples.GE2016only.castsOnly.mat
% GE2016 = samplesGE2016castsOnly;

% load gcms_greenedge_proc
% DATA = [data DMSOUT];
load ALLSURF

%% GENERAL SETTINGS
h = 'z10'; % hBD, mld, isolu, z60, z10
varlist = {'dms' 'fdmsW97c24' 'OWD'...
    sprintf('idms_%s',h) sprintf('idmspt_%s',h) sprintf('ddratio_%s',h)...
    sprintf('icp_%s',h) sprintf('iTchla_%s',h) sprintf('ccratio_%s',h)};
titles = {'DMS (nM)' 'FDMS (µmol m^{-2} d^{-1})'...
    'Open Water Days'...
    '\SigmaDMS (µmol m^{-2})' '\SigmaDMSPt (µmol m^{-2})' '\SigmaDMS:\SigmaDMSPt'...
    '\Sigmac_{p} (m^{-1} m)' '\SigmaTChla (mg m^{-2})' '\Sigmac_{p}:\SigmaTChla'};
icol.fdms = strcmp(headALLSURF,'fdmsW97c24');
icol = strcmp(headALLSURF,'dms');
icolpt = strcmp(headALLSURF,'dmspt');
% colors = brewermap(11,'Blues');
% ccmat = colors([6 8 11],:);
% colors = brewermap(11,'Spectral');
% ccmat = colors([10 9 3],:);
% ccmat = colors([11 9 2],:);
colors = brewermap(21,'Spectral');
ccmat = colors([20 17 6],:);
xtl = {'ICE' 'MIZ' 'OW'};
ncols = 3;
nrows = 3;
x0offset = 1.5;
y0offset = 1.5;
xtoffset = 1;
ytoffset = 1;
xspace = 1.3;
yspace = 1.1;
pwidth = 3;
pheight = 4;
twidth = ncols*(pwidth+xspace)+x0offset;
theight = nrows*(pheight+yspace)+y0offset;
ymintomax = 0.15;
if strcmp(h,'hBD')
    ymaxs = [20 20 25 250 2500 0.33 10 30  0.9];
elseif strcmp(h,'mld')
    ymaxs = [20 20 25 250 2500 0.33 10 30  0.9];
elseif strcmp(h,'z10')
    ymaxs = [20 20 23 200 2000 0.4 8 20  0.8];
%     yletters = [15 15 11 150 1500 0.3 6 15 0.6];
elseif strcmp(h,'isolu')
    ymaxs = [20 20 25 900 9000 0.33 30 100  0.9];
elseif strcmp(h,'z60')
    ymaxs = [20 20 25 900 9000 0.33 30 100  0.9];
end
yletters = [17.5 17.5 17.5 175 1750 0.35 7 17 0.7];
pletters = {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i'};
ymins = -ymaxs*ymintomax;
ymins(3) = -25; % Specific for OWD
ytexts = ymins + 0.07*(ymaxs-ymins);

%% Classify by ice cover
icemean = nanmean(ALLSURF(:,11:13),2); % 3 day ice mean
icemin = nanmin(ALLSURF(:,11:13),[],2);
icemax = nanmax(ALLSURF(:,11:13),[],2);
icecrit1 = 0.10;
icecrit2 = 0.85;
icebin = 2*ones(size(icemax));
icebin(icemax<icecrit1) = 3; % Open: Persistent <10% ice
icebin(icemin>icecrit2) = 1; % Covered: Persistent >80% ice (85% gives similar results)

% Set to NaN OWD data where no FDMS is available
ALLSURF(isnan(ALLSURF(:,32)),strcmp(headALLSURF,'OWD')) = nan;

% % Set to NaN icp and iTChla data where no dmspt is available
% ALLSURF(isnan(ALLSURF(:,31)),45:49) = nan;
% ALLSURF(isnan(ALLSURF(:,31)),51:5) = nan;

% OR: Set to NaN icp and iTChla data for leg 1a
ALLSURF(ALLSURF(:,1)<400,45:49) = nan;
ALLSURF(ALLSURF(:,1)<400,51:55) = nan;

%% Append ratios
ddratio = ALLSURF(:,strcmp(headALLSURF,sprintf('idms_%s',h)))./ALLSURF(:,strcmp(headALLSURF,sprintf('idmspt_%s',h)));
ccratio = ALLSURF(:,strcmp(headALLSURF,sprintf('icp_%s',h)))./ALLSURF(:,strcmp(headALLSURF,sprintf('iTchla_%s',h)));
ALLSURF = [ALLSURF ddratio ccratio];
headALLSURF = [headALLSURF; sprintf('ddratio_%s',h); sprintf('ccratio_%s',h)];

%% Correct DMS flux multiplying by SIC
dayice = nanmean([ALLSURF(:,strcmp(headALLSURF,'SICday')) ALLSURF(:,strcmp(headALLSURF,'IceConcentration_percent'))/100],2);
Fice = ALLSURF(:,strcmp(headALLSURF,'fdmsW97c24'));
nanmean(Fice)
% ALLSURF(:,strcmp(headALLSURF,'fdmsW97c24')) = Fice.*(1-dayice);
ALLSURF(:,strcmp(headALLSURF,'fdmsW97c24')) = Fice.*(1-icemean);
nanmean(ALLSURF(:,strcmp(headALLSURF,'fdmsW97c24')))

%% Boxplots by ice cover

figure(33), clf
set(gcf,'units','centimeters','position',[1 2 twidth theight])
% set(gcf,'units','centimeters','position',[1 2 15 20])

for j = 1:nrows
    for k = 1:ncols
        jk = (k-1)*nrows+j;
        sposition = [x0offset+(j-1)*(pwidth+xspace) theight-k*(pheight+yspace) pwidth pheight];
        sposition([1 3]) = sposition([1 3])/twidth;
        sposition([2 4]) = sposition([2 4])/theight;
        subplot('position',sposition);
        %         subplot(nrows,ncols,jk)
        icol = strcmp(headALLSURF,varlist{jk});
        xplot = ALLSURF(:,icol);
        % Ice cover correction for FDMS only
        if strcmp(varlist{jk},'fdmsW97c24') || strcmp(varlist{jk},'fdmsW97c24')
            xplot = xplot.*(1-icemean);
        end
        % Traditional
        boxplot(xplot,icebin,'colors',ccmat,'symbol','.w','outliersize',0.1,'notch','on',...
            'boxstyle','outline','widths',.5,'plotstyle','traditional','medianstyle','line'), hold on
        %         % Compact
        %         boxplot(xplot,icebin,'colors',ccmat,'symbol','+k',...
        %             'plotstyle','compact','width',.5), hold on
        for ii = 1:3
            plot(ii+[-.15 .15],[1 1]*nanmedian(xplot(icebin==ii)),'-','markersize',10,...
                'color',ccmat(ii,:),'linewidth',2)
            plot(ii,nanmean(xplot(icebin==ii)),'o','markersize',4,...
                'markeredgecolor',ccmat(ii,:),'markerfacecolor','w','linewidth',1.5)
            text(ii-0.1,ytexts(jk),sprintf('%i',sum(~isnan(xplot(icebin==ii)))),'fontsize',5,'color',[.5 .5 .5])
            % text(ii-0.2,ytexts(jk),sprintf('N=%i',sum(~isnan(xplot(icebin==ii)))),'fontsize',5,'color',[.5 .5 .5])
        end
        title(titles(jk),'fontsize',9,'fontweight','bold')
        text(0.7,yletters(jk),sprintf('%s',pletters{jk}),'fontsize',9,'color','k','fontweight','bold')
        ylim([ymins(jk) ymaxs(jk)])
        set(gca,'tickdir','in','ticklength',[.02 .02],'xtick',1:max(icebin),...
            'xticklabel',{''},'fontsize',8,'linewidth',0.75,...
            'xlim',[0.5 3.5])
        if k == 3
            set(gca,'xticklabel',xtl)
        end
    end
end

% NOTE
% boxstyle filled vs outline (allows for width to be set)
% plotstyle compact vs traditional
% symbol refers only to outliers


%% Plot without categorization, x = OWD

% % Sort by lon, lat and doy in this order
% [GE2016, indexs] = sortrows(GE2016,[7 6]);
% ALLSURF = ALLSURF(indexs,:);

%%
x = ALLSURF(:,9); sum(~isnan(x)) % OWD
% sic = nanmean([ALLSURF(:,4) 100*ALLSURF(:,13)],2); sum(~isnan(sic)) 
sic = nanmean([100*ALLSURF(:,11:13)],2); sum(~isnan(sic)) 

%% Interpolations of owd and sic: interp 2D and 3D a bit complicate.
% Try instead with nearest neighbor: not good idea

% lon = GE2016(:,7);
% lat = GE2016(:,6);
% doy = GE2016(:,2);
% dlon = 0.1;
% dlat = 0.05;
% ddoy = 1;
% 
% %% indices of x points to be "interpolated"
% ix = find(isnan(x))';
% for j = ix
%     ilon = abs(lon(j) - lon) <= dlon;
%     ilat = abs(lat(j) - lat) <= dlat;
%     idoy = abs(doy(j) - doy) <= ddoy;
%     ij = ilon & ilat & idoy;
%     sum(ij)
%     if sum(ij) ~= 0
%         x(j) = nanmean(x(ij));
%     end
% end
% sum(~isnan(x))
% 
% %% indices of owd points to be "interpolated"
% is = find(isnan(sic))';
% for j = is
%     ilon = abs(lon(j) - lon) <= dlon;
%     ilat = abs(lat(j) - lat) <= dlat;
%     idoy = abs(doy(j) - doy) <= ddoy;
%     ij = ilon & ilat & idoy;
%     sum(ij)
%     if sum(ij) ~= 0
%         sic(j) = nanmean(x(ij));
%     end
% end
% sum(~isnan(sic))

%%
figure(333), clf

subplot(2,1,1)
hold on
plot(x, sic,'ob')
plot(x, ALLSURF(:,17)*10+10,'+r')
plot(x, ALLSURF(:,18),'+c')
plot(x, ALLSURF(:,14)*10,'xk')
plot(x, ALLSURF(:,10),'.y','markersize',20)
plot(x, ALLSURF(:,32),'.m','markersize',20)
% plot(x, ALLSURF(:,8)*100,'*g','markersize',10) % Atl vs Arc clustering not too clear

subplot(2,1,2)
hold on
plot(x, -ALLSURF(:,6),'.k','markersize',20)
plot(x, -ALLSURF(:,24),'.k','markersize',10)
plot(x, -ALLSURF(:,30),'.g','markersize',20)
plot(x, -ALLSURF(:,3),'.y','markersize',10)

%% Plot without categorization, x = 3-day ice cover. Not too clear

% x = nanmean(ALLSURF(:,11:13),2); % OWD
% % sic = 100*nanmean([ALLSURF(:,4) 100*ALLSURF(:,13)],2);
% sic = nanmean([100*ALLSURF(:,11:13)],2);
% 
% figure(334), clf
% hold on
% plot(x, ALLSURF(:,9),'ob')
% plot(x, ALLSURF(:,17)*10,'+r')
% plot(x, ALLSURF(:,18),'+c')
% plot(x, ALLSURF(:,14),'xk')
% plot(x, ALLSURF(:,10),'.y','markersize',20)
% plot(x, -ALLSURF(:,6),'.k','markersize',20)
