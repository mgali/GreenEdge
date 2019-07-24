clear, clc

% Variables included in analysis
headX = {'stn' 'depth' 'temp' 'sal' 's-t' 'O2' 'cp' 'DMS'};

%% Load ATOS07 data

load('profile_ATOS07.mat')
D.A07 = DATA_CTD(~isnan(DATA_CTD(:,8)),:);

% remove Greenland current stations
% D.A07(D.A07(:,1)<9,:) = nan;
% D.A07(D.A07(:,1)>=9,:) = nan;

%% Load and select GE data
% Need to convert O2 from mg/L to umol/kg

load('surf_profile_ALL_4plot.mat','SURF','headSURF','DATA_CTD','headDATA_CTD')
iselect = [6 13 48 49 50 52 56 45]; % stn z temp sal sigmat O2 cpsmooth1 DMS
D.GE = DATA_CTD(:,iselect);

% remove Greenland current stations
% D.GE(D.GE(:,1)<500,:) = nan;
% D.GE(D.GE(:,1)>=500,:) = nan;

%% DMS vs cp

% figure(1), clf, hold on
% x1 = real(log10(D.A07(:,7)));
% y1 = real(log10(D.A07(:,8)));
% p1 = plot(x1,y1,'o','markersize',12,'markeredgecolor','k','markerfacecolor',[1 .7 .3]);
% x2 = real(log10(D.GE(:,7)));
% y2 = real(log10(D.GE(:,8)));
% p2 = plot(x2,y2,'o','markersize',12,'markeredgecolor','k','markerfacecolor',[0 .7 .3]);
%
% x = [x1; x2];
% y = [y1; y2];
% excl = isnan(x) | isnan(y) | abs(x)==Inf;
% x(excl) = [];
% y(excl) = [];
%
% xlabel('c_{p} (m^{-1})','fontsize',16)
% ylabel('DMS (nM)','fontsize',16)
% set(gca,'tickdir','out','fontsize',14,...
%     'xticklabel',[.003 .01 .03 .1 .3 1 3],...
%     'yticklabel',[.03 .1 .3 1 3 10 30 100])
% text(-0.5,-0.7,sprintf('r_{lin} = %0.2f',corr(x,y,'rows','pairwise')),'fontsize',16)
% text(-0.5,-1,sprintf('r_{log} = %0.2f',corr(x,y,'rows','pairwise')),'fontsize',16)


% Fit

% xy = sortrows(10.^([x y]),1);
% x = xy(:,1); y = xy(:,2);
% mdl = fit(x,y,'poly1');
% ci.obs = predint(mdl,x,0.68,'observation','off');
% yp = x*mdl.p1 + mdl.p2;
%
% plot(log10(x),log10(yp),'-k','linewidth',2)
% ci.obs(ci.obs<0) = nan;
% plot(log10(x),log10(ci.obs),'-k')


%% DMS vs cp by transects (GE) and water masses (A07)

% Create groups
D.A07(D.A07(:,1)<9,1) = 1; % East Greenland Current
D.A07(D.A07(:,1)==9|D.A07(:,1)==12|D.A07(:,1)==20|D.A07(:,1)==23|D.A07(:,1)==27|D.A07(:,1)==33,1) = 2; % AW
D.A07(D.A07(:,1)==18|D.A07(:,1)==19|D.A07(:,1)==26|D.A07(:,1)==30|D.A07(:,1)>=36,1) = 3; % PSW
D.GE(:,1) = floor(D.GE(:,1)/100);
D.TOT = [D.A07; D.GE];
stn = D.TOT(:,1);

clist = flip(brewermap(8,'Paired'));
clist(4,:) = [];
symbols = {'o' 's' 'o' 's' 'o' 's' 'o' 's' 'o'};

figure(2), clf, hold on
for j = 1:7
    TMP = D.TOT(stn==j,:);
    x1 = real(log10(TMP(:,7)));
    y1 = real(log10(TMP(:,8)));
    plot(x1,y1,symbols{j},'markersize',10,'markeredgecolor','k','markerfacecolor',clist(j,:));
end
legend({'A-EGC' 'A-AW' 'A-PSW' 'G400' 'G500' 'G600' 'G700'},'location','best')

% x = [x1; x2];
% y = [y1; y2];
% excl = isnan(x) | isnan(y) | abs(x)==Inf;
% x(excl) = [];
% y(excl) = [];

xlabel('c_{p} (m^{-1})','fontsize',16)
ylabel('DMS (nM)','fontsize',16)
set(gca,'tickdir','out','fontsize',14,...
    'xticklabel',[.003 .01 .03 .1 .3 1 3],...
    'yticklabel',[.03 .1 .3 1 3 10 30 100])
% text(-0.5,-0.7,sprintf('r_{lin} = %0.2f',corr(x,y,'rows','pairwise')),'fontsize',16)
% text(-0.5,-1,sprintf('r_{log} = %0.2f',corr(x,y,'rows','pairwise')),'fontsize',16)

