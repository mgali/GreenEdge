clear, clc

load ALLSURF
load samples.GE2016only.castsOnly.mat
GE2016 = samplesGE2016castsOnly;

%% --------------- Check consistency in ALLSURF matrix ---------------

%% Achim's vs. my MLD
figure(101), clf
% x = ALLSURF(:,strcmp(headALLSURF,'mld125'));
x = ALLSURF(:,strcmp(headALLSURF,'mld03'));
y = ALLSURF(:,strcmp(headALLSURF,'MixedLayerDepth_0.1gkg-1-criterion_m'));
scatter(x,y,80), hold on
plot([0 80],[0 80])
% xlabel('MLD125')
xlabel('MLD03')
ylabel('MLD100')

%% Achim's OWD vs. 3-day SIC
figure(102), clf
x = ALLSURF(:,strcmp(headALLSURF,'IceConcentration_percent'));
y = nanmean(ALLSURF(:,strcmp(headALLSURF,'SICday')),2);
scatter(x,y,80)
xlabel('SIC % Achim')
ylabel('SIC day')

%% Achim's OWD vs. 3-day SIC
% figure(103), clf
% x = ALLSURF(:,strcmp(headALLSURF,'OWD'));
% y = nanmean(ALLSURF(:,strcmp(headALLSURF,'SICm2d')|strcmp(headALLSURF,'SICm1d')|strcmp(headALLSURF,'SICday')),2);
% scatter(x,y,80)
% xlabel('OWD')
% ylabel('SIC 3 days')


%% --------------- Explore relationships in ALLSURF matrix ---------------

%% Integrated stocks
% layer = 'mld';
% layer = 'hBD';
% layer = 'isolu';
layer = 'z60';
figure(104), clf
x = ALLSURF(:,strcmp(headALLSURF,strcat('icp','_',layer)));
y = nanmean(ALLSURF(:,strcmp(headALLSURF,strcat('idmspt','_',layer))),2);
plot(x,y/10,'.r','markersize',40), hold on
y = nanmean(ALLSURF(:,strcmp(headALLSURF,strcat('idms','_',layer))),2);
plot(x,y,'ob','markersize',10)
xlabel('Cp integral, n.d.','fontsize',16)
ylabel('DMS or DMSPt/10 integral, µmol m^{-2}','fontsize',16)
title(sprintf('Integrated layer: %s',layer),'fontsize',16)

%% Integrated stocks map
cmap = parula;
% layer = 'mld';
% layer = 'hBD';
layer = 'isolu';
% layer = 'z60';
intvar = 'cp';
figure(104), clf
st = ALLSURF(:,strcmp(headALLSURF,'station'));
lat = GE2016(:,6);
lon = GE2016(:,7);
z = ALLSURF(:,strcmp(headALLSURF,strcat('i',intvar,'_',layer)));
scatter(lon,lat,80,z,'filled')
xlim([-65 -56])
ylim([67.5 71])
xlabel('Longitude E','fontsize',16)
ylabel('Latitude N','fontsize',16)
title(sprintf('Integrated layer: %s',layer),'fontsize',16)
colormap(cmap)
c = colorbar;
ylabel(c,sprintf('%s integral, n.d.',intvar),'fontsize',16)
set(c,'fontsize',14)

