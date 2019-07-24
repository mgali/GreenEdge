% PLOT HISTOGRAMS OF IN SITU PAR TIME SERIES
clc
clear


%% Ice camp 2015
load('~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2015-ICECAMP_cyber_fractSZA41_44.mat')
% load('~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2015-ICECAMP_cyber_fractSZA42_45.mat')
% load('~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2015-ICECAMP_cyber_fractSZA40_45.mat')
% load('~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2015-ICECAMP_cyber_fract43.mat')
% load('~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2015-ICECAMP_cyber_fract42.mat')
par1 = par_ts.data;

% Record xlimits
dv = datevec(par_ts.mtimeUTC);
doy1 = yearday(dv(:,3),dv(:,2),dv(:,1),dv(:,4),dv(:,5),dv(:,6));
xlim1 = [nanmin(doy1); nanmax(doy1)];

% SBDART DATA
% load('~/Desktop/GreenEdge/Irradiance/SBDART_LUTs_outputs/PAR_SBDART_GE2015-ICECAMP_SurfAlb_v1_1h_fakeIce.mat')
% load('~/Desktop/GreenEdge/Irradiance/SBDART_LUTs_outputs/PAR_SBDART_GE2015-ICECAMP_SurfAlb_v1_1h_realIce.mat')
load('~/Desktop/GreenEdge/Irradiance/SBDART_LUTs_outputs/PAR_SBDART_GE2015-ICECAMP_SurfAlb_v1_1h_consensus.mat')
dv = datevec(par_ts.mtimeUTC);
doy11 = yearday(dv(:,3),dv(:,2),dv(:,1),dv(:,4),dv(:,5),dv(:,6));
par11 = par_ts.data;
uvr11 = par_ts.dataUVR;

%% Ice camp 2016
load('~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2016-ICECAMP_cyber.mat')
par2 = par_ts.data;

% Record xlimits
dv = datevec(par_ts.mtimeUTC);
doy2 = yearday(dv(:,3),dv(:,2),dv(:,1),dv(:,4),dv(:,5),dv(:,6));
xlim2 = [nanmin(doy2); nanmax(doy2)];

% SBDART DATA
% load('~/Desktop/GreenEdge/Irradiance/SBDART_LUTs_outputs/PAR_SBDART_GE2016-ICECAMP_SurfAlb_v1_1h_fakeIce.mat')
% load('~/Desktop/GreenEdge/Irradiance/SBDART_LUTs_outputs/PAR_SBDART_GE2016-ICECAMP_SurfAlb_v1_1h_realIce.mat')
load('~/Desktop/GreenEdge/Irradiance/SBDART_LUTs_outputs/PAR_SBDART_GE2016-ICECAMP_SurfAlb_v1_1h_consensus.mat')
dv = datevec(par_ts.mtimeUTC);
doy22 = yearday(dv(:,3),dv(:,2),dv(:,1),dv(:,4),dv(:,5),dv(:,6));
par22 = par_ts.data;
uvr22 = par_ts.dataUVR;

%% Common x limits and x annotation

xmin = floor(nanmin([xlim1; xlim2]));
xmax = ceil(nanmax([xlim1; xlim2]));
xax = xmin:xmax;
xt = xax(~mod(xax,5));

%%
figure(545), clf
set(gcf,'units','centimeters','position',[5 10 30 14.5])


subplot(2,1,1), hold on
plot(doy11, par11, '-r', 'linewidth', 0.5), hold on
plot(doy1, par1, '.k', 'markersize', 0.1)
set(gca,'units','centimeters','position',[2 8 27 6],...
    'tickdir','out','xtick',xt,'xticklabel','')
ylabel('PAR (µmol photons m^{-2} s^{-1})')
xlim([xmin xmax])
ylim([0 2400])
text(180,2000,'- In situ','color','k')
text(180,1800,'- SBDART/Aqua','color','r')

subplot(2,1,2), hold on
plot(doy22, par22, '-r', 'linewidth', 0.5), hold on
plot(doy2, par2, '-k', 'linewidth', 0.5)
set(gca,'units','centimeters','position',[2 2 27 5],...
    'tickdir','out','xtick',xt)
xlabel('DOY')
ylabel('PAR (µmol photons m^{-2} s^{-1})')
xlim([xmin xmax])
ylim([0 2000])

