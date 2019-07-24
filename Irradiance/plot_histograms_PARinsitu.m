% PLOT HISTOGRAMS OF IN SITU PAR TIME SERIES
clc
clear

figure(454), clf
set(gcf,'units','centimeters','position',[5 10 30 10])

% General definitions
bs = 100;
bins = (0:bs:2500)';
colors = {'r' [0 .6 0] 'b'};
leg = ({'Ice camp 2015' 'Ice camp 2016' 'Amundsen 2016'})';
Y = nan(length(bins),3);

%%

% Ice camp 2015
load('~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2015-ICECAMP_cyber_fractSZA41_44.mat')
% load('~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2015-ICECAMP_cyber_fractSZA42_45.mat')
% load('~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2015-ICECAMP_cyber_fractSZA40_45.mat')
% load('~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2015-ICECAMP_cyber_fract43.mat')
% load('~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2015-ICECAMP_cyber_fract42.mat')
x = par_ts.data;
Y(:,1) = histc(x,bins);

% Ice camp 2016
load('~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2016-ICECAMP_cyber.mat')
x = par_ts.data;
Y(:,2) = histc(x,bins);

% Amundsen 2016
load('~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2016-Amundsen_cyber.mat')
x = par_ts.data;
cutx = datenum([2016 6 9 0 0 0]);
x(par_ts.mtimeUTC < cutx) = nan;
Y(:,3) = histc(x,bins);

% Normalize to their counts
ycounts = sum(Y);
Yn = Y./(ones(size(Y,1),1)*ycounts);

for j = 1:size(Yn,2)
    leg{j} = sprintf('%s, n = %i', leg{j}, ycounts(j));
end

%%
subplot(1,2,1), hold on
set(gca,'units','centimeters','position',[2 1.5 13 8])
for j = 1:size(Yn,2)
    plot(bins+bs/2,Yn(:,j),'-+','color',colors{j},'linewidth',1)
end

set(gca,'xtick',bins(~mod(bins,200)),'xticklabel',bins(~mod(bins,200)),...
    'tickdir','out','xscale','lin','yscale','lin')
xlabel('PAR (µmol photons m^{-2} s^{-1})')
ylabel('Normalized frequency')
grid on
xlim([nanmin(bins) nanmax(bins)])
legend(leg,'location','best')
legend boxoff

%%
subplot(1,2,2), hold on
set(gca,'units','centimeters','position',[16.5 1.5 13 8])
for j = 1:size(Yn,2)
    plot(bins+bs/2,Yn(:,j),'-+','color',colors{j},'linewidth',1)
    leg{j} = sprintf('%s, n = %i', leg{j}, ycounts(j));
end

set(gca,'xtick',bins(~mod(bins,200)),'xticklabel',bins(~mod(bins,200)),...
    'tickdir','out','xscale','lin','yscale','log')
xlabel('PAR (µmol photons m^{-2} s^{-1})')
% ylabel('Frequency')
xlim([nanmin(bins) nanmax(bins)])
grid on
% legend(leg,'location','best')
legend boxoff



