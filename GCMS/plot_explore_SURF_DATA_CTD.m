clear, clc

load('surf_profile_ALL_4plot.mat','SURF','headSURF','DATA_CTD','headDATA_CTD')


%% --------------- Explore relationships in DATA_CTD matrix ---------------

%% DMS vs X all data. X is cp or NO3

% figure(1), clf
% x = log10(DATA_CTD(:,strcmp(headDATA_CTD,'cpsmooth1'))); % cpsmooth1, NO3
% y = log10(DATA_CTD(:,strcmp(headDATA_CTD,'dms_consens_cf68')));
% z = DATA_CTD(:,13);
% s = scatter(x,y,80,-z,'filled');
% xlabel('c_{p} (m^{-1})','fontsize',16)
% % xlabel('NO_{3} (µM)','fontsize',16)
% ylabel('[DMS] (nM)','fontsize',16)
% set(gca,'tickdir','out','fontsize',14,...
%     'xticklabel',[.003 .01 .03 .1 .3 1 3],... % 100* if nitrate
%     'yticklabel',[.03 .1 .3 1 3 10 30 100])
% text(-0.3,-0.5,sprintf('r_{lin} = %0.2f',corr(10.^x,10.^y,'rows','pairwise')),'fontsize',16)
% text(-0.3,-1,sprintf('r_{log} = %0.2f',corr(x,y,'rows','pairwise')),'fontsize',16)
% 
% c = colorbar;
% cmap = flip(parula);
% colormap(cmap)
% ylabel(c,'Depth (m)','fontsize',16)
% set(c,'fontsize',14)

%% DMS vs CDOM all data

% figure(1), clf
% x = log10(DATA_CTD(:,strcmp(headDATA_CTD,'CDOM')));
% y = log10(DATA_CTD(:,strcmp(headDATA_CTD,'dms_consens_cf68')));
% z = DATA_CTD(:,13);
% s = scatter(x,y,80,-z,'filled');
% xlabel('CDOM (units?)','fontsize',16)
% ylabel('[DMS] (nM)','fontsize',16)
% set(gca,'tickdir','out','fontsize',14,...
%     'xticklabel',0.2:0.05:0.5,...
%     'yticklabel',[.03 .1 .3 1 3 10 30 100]) % 
% text(0.25,-0.5,sprintf('r_{lin} = %0.2f',corr(10.^x,10.^y,'rows','pairwise')),'fontsize',16)
% text(0.25,-1,sprintf('r_{log} = %0.2f',corr(x,y,'rows','pairwise')),'fontsize',16)
% 
% c = colorbar;
% cmap = flip(parula);
% colormap(cmap)
% ylabel(c,'Depth (m)','fontsize',16)
% set(c,'fontsize',14)

%% DMS vs O2 all data LOG

% figure(1), clf
% x = log10(DATA_CTD(:,strcmp(headDATA_CTD,'O2')));
% y = log10(DATA_CTD(:,strcmp(headDATA_CTD,'dms_consens_cf68')));
% z = DATA_CTD(:,13);
% s = scatter(x,y,80,-z,'filled');
% xlabel('O_{2} (mg L^{-1})','fontsize',16)
% ylabel('[DMS] (nM)','fontsize',16)
% set(gca,'tickdir','out','fontsize',14,...
%     'xticklabel',[6.3 7.1 7.9 8.9 10],...
%     'yticklabel',[.03 .1 .3 1 3 10 30 100]) % 
% text(0.95,-0.5,sprintf('r_{lin} = %0.2f',corr(10.^x,10.^y,'rows','pairwise')),'fontsize',16)
% text(0.95,-1,sprintf('r_{log} = %0.2f',corr((x),(y),'rows','pairwise')),'fontsize',16)
% 
% c = colorbar;
% cmap = flip(parula);
% colormap(cmap)
% ylabel(c,'Depth (m)','fontsize',16)
% set(c,'fontsize',14)

%% DMS vs O2 all data LIN

% figure(1), clf
% x = DATA_CTD(:,strcmp(headDATA_CTD,'O2'));
% y = DATA_CTD(:,strcmp(headDATA_CTD,'dms_consens_cf68'));
% z = DATA_CTD(:,13);
% s = scatter(x,y,80,-z,'filled');
% xlabel('O_{2} (mg L^{-1})','fontsize',16)
% ylabel('[DMS] (nM)','fontsize',16)
% set(gca,'tickdir','out','fontsize',14,...
%     'xticklabel',6.5:0.5:10,...
%     'yticklabel',0:10:80)
% text(9,70,sprintf('r_{lin} = %0.2f',corr(x,y,'rows','pairwise')),'fontsize',16)
% text(9,60,sprintf('r_{log} = %0.2f',corr(log10(x),log10(y),'rows','pairwise')),'fontsize',16)
% 
% c = colorbar;
% cmap = flip(parula);
% colormap(cmap)
% ylabel(c,'Depth (m)','fontsize',16)
% set(c,'fontsize',14)


%% Multivariate

% % NO VARIABLE WORTH ADDING AFTER CP

yr = DATA_CTD(:,strcmp(headDATA_CTD,'dms_consens_cf68'));
% istepw = [13 48:54 56];
istepw = [52:54 56];
Xr = DATA_CTD(:,istepw);
headSTEPW = headDATA_CTD(istepw);
% [b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(Xr,yr,'penter',0.05,'premove',0.10);
[b,bi,r,ri,st] = regress(yr,Xr);

%% --------------- Explore relationships in SURF matrix ---------------

%% Classify by ice cover (open receding covered)

icemean = nanmean(SURF(:,8:10),2);
icemin = nanmin(SURF(:,8:10),[],2);
icemax = nanmax(SURF(:,8:10),[],2);
icecrit1 = 0.10;
icecrit2 = 0.80;
icebin = 2*ones(size(icemax));
icebin(icemax<icecrit1) = 1; % Open: Persistent <10% ice
icebin(icemin>icecrit2) = 3; % Covered: Persistent >80% ice
corrflux = SURF(:,16).*(1-icemean);


%% DMS vs X surface

% figure(2), clf
% x = log10(DATA_CTD(DATA_CTD(:,13)<15,end));
% y = log10(DATA_CTD(DATA_CTD(:,13)<15,idms));
% z = DATA_CTD(DATA_CTD(:,13)<15,13);
% s = scatter(x,y,80,-z,'filled');
% xlabel('log10(c_{p})','fontsize',16)
% ylabel('log10(DMS)','fontsize',16)
% % set(gca,'tickdir','out','fontsize',14,...
% %     'xticklabel',[.003 .01 .03 .1 .3 1 3],...
% %     'yticklabel',[.03 .1 .3 1 3 10 30 100])
% % text(-0.5,-0.5,sprintf('r = %0.2f',corr(x,y,'rows','pairwise')),'fontsize',16)
% % text(-0.5,-1,sprintf('r_{S} = %0.2f',corr(x,y,'rows','pairwise')),'fontsize',16)
% 
% c = colorbar;
% cmap = flip(parula);
% colormap(cmap)
% ylabel(c,'Depth (m)','fontsize',16)
% set(c,'fontsize',14)

%% DMS vs MLD03, by SIC

% cmap = brewermap(11,'Blues');
% 
% figure(101), clf
% x = SURF(:,strcmp(headSURF,'mld03'));
% y = SURF(:,strcmp(headSURF,'dms'));
% z = SURF(:,strcmp(headSURF,'SICday'));
% scatter(x,y,80,z,'filled','markeredgecolor','k')
% colormap(cmap)
% xlabel('MLD03')
% ylabel('DMS')
% c = colorbar;
% ylabel(c,'SIC (same day)','fontsize',16)
% set(c,'fontsize',14)

%% DMS vs MLD125, by SIC

% cmap = brewermap(11,'Blues');
% 
% figure(102), clf
% x = SURF(:,strcmp(headSURF,'mld125'));
% y = SURF(:,strcmp(headSURF,'dms'));
% z = SURF(:,strcmp(headSURF,'SICday'));
% scatter(x,y,80,z,'filled','markeredgecolor','k')
% colormap(cmap)
% xlabel('MLD125')
% ylabel('DMS')
% c = colorbar;
% ylabel(c,'SIC (same day)','fontsize',16)
% set(c,'fontsize',14)

%% DMS vs cp, by X (WS, SIC)

% cmap = brewermap(11,'Reds');
cmap = flip(brewermap(11,'Blues'));

figure(102), clf
x = SURF(:,strcmp(headSURF,'cpsmooth1'));
y = SURF(:,strcmp(headSURF,'dms'));
z = SURF(:,strcmp(headSURF,'SICday'));
scatter(x,y,80,z,'filled','markeredgecolor','k')
colormap(cmap)

xlabel('c_{p} (m^{-1})','fontsize',14)
ylabel('DMS (nM)','fontsize',14)
c = colorbar;
% ylabel(c,'Wind speed (previous 72h)','fontsize',16)
ylabel(c,'SIC (same day)','fontsize',16)
set(c,'fontsize',14)

%% DMS vs MLD125, by cp

% figure(103), clf
% x = SURF(:,strcmp(headSURF,'mld125'));
% y = SURF(:,strcmp(headSURF,'dms'));
% z = SURF(:,strcmp(headSURF,'cp'));
% scatter(x,y,80,z,'filled')
% xlabel('MLD125')
% ylabel('DMS')
% c = colorbar;
% ylabel(c,'cp (m^{-1})','fontsize',16)
% set(c,'fontsize',14)

%% DMS vs N2_MLD03, by SIC

% figure(104), clf
% x = SURF(:,strcmp(headSURF,'N2max03'));
% y = SURF(:,strcmp(headSURF,'dms'));
% z = SURF(:,strcmp(headSURF,'SICday'));
% scatter(x,y,80,z,'filled')
% xlabel('N2max-MLD03')
% ylabel('DMS')
% c = colorbar;
% colormap parula
% ylabel(c,'SIC (same day)','fontsize',16)
% set(c,'fontsize',14)

%% DMS vs N2_MLD125, by SIC

% figure(105), clf
% x = SURF(:,strcmp(headSURF,'N2max125'));
% y = SURF(:,strcmp(headSURF,'dms'));
% z = SURF(:,strcmp(headSURF,'SICday'));
% scatter(x,y,80,z,'filled')
% xlabel('N2max-MLD125')
% ylabel('DMS')
% c = colorbar;
% colormap parula
% ylabel(c,'SIC (same day)','fontsize',16)
% set(c,'fontsize',14)

%% MLD03 vs N2_MLD03, by SIC

% figure(106), clf
% x = SURF(:,22);%strcmp(headSURF,'N2max03'));
% y = SURF(:,strcmp(headSURF,'mld03'));
% z = SURF(:,strcmp(headSURF,'SICday'));
% scatter(x,y,80,z,'filled')
% xlabel('N2max-MLD03')
% ylabel('MLD03')
% c = colorbar;
% colormap parula
% ylabel(c,'SIC (same day)','fontsize',16)
% set(c,'fontsize',14)

%% MLD125 vs N2_MLD125, by SIC

% figure(107), clf
% x = SURF(:,23);%strcmp(headSURF,'N2max125'));
% y = SURF(:,strcmp(headSURF,'mld125'));
% z = SURF(:,strcmp(headSURF,'SICday'));
% scatter(x,y,80,z,'filled')
% xlabel('N2max-MLD125')
% ylabel('MLD125')
% c = colorbar;
% colormap parula
% ylabel(c,'SIC (same day)','fontsize',16)
% set(c,'fontsize',14)


%% Multivariate

% % SELECTED VARIABLES ARE COINCIDENTAL RATHER THAN DRIVING FACTORS
% % For example, WS and SST are higher in ice-free zone where surface DMS
% % is higher. But not because of higher SST or higher WS.

% yr = SURF(:,strcmp(headSURF,'dms'));
% istepw = [10 11 14 15 18:23 25 26];
% Xr = SURF(:,istepw);
% headSTEPW = headSURF(istepw);
% [b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(Xr,yr,'penter',0.05,'premove',0.10);
% % [b,bi,r,ri,st] = regress(yr,Xr);
% 
% ifinal = istepw(inmodel);
% vars_in = headSURF(ifinal)
% b_in = b(inmodel)
% se_in = se(inmodel)

