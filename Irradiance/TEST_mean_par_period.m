% Test mean_par_period.m on 2016 ice camp PAR data, plot
clc, clear

%% CASE 1, ice camp 2016

load ~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2016-ICECAMP_cyber.mat
sample_list = dlmread('input.noclim.IceCamp20152016.txt');

% General definitions
tsearch_h = 3.1;
pXh = 48;
UTC_diff = 4;
lon = -63.78953;
real_diff = -lon/15;
noonLOC = 0.5;
real_noonLOC = noonLOC + (real_diff-UTC_diff)/24;
tsamLOC = 10/24;

% Preallocate and fill
PAR0PLUS = nan(size(sample_list,1),19);

for id = 1:size(sample_list,1)
%     if mod(id,2)
%         tsamplUTC = datenum(sample_list(id,5),sample_list(id,3),sample_list(id,4),14,0.5,0);
%     else
%         tsamplUTC = datenum(sample_list(id,5),sample_list(id,3),sample_list(id,4),18,0.5,0);
%     end
    tsamplUTC = datenum(sample_list(id,5),sample_list(id,3),sample_list(id,4),14,0.5,0);
    par0plus = mean_par_period(par_ts, tsamplUTC, tsearch_h, pXh, lon, UTC_diff, 1);
    PAR0PLUS(id,:) = (cell2mat(struct2cell(par0plus)))';
end

% --------- Plot original time series along with extracted values ---------

PAR2PLOT = PAR0PLUS(:,[1 2 7 8 9 10 12 13 14]);

% Select time axis min and max
tmin_plot = yearday(1,5,2016,0,0,0);
tmax_plot = yearday(11,5,2016,0,0,0);

% Define various time axes
tref = datenum([2016 1 1 0 0 0]);
tnew = par_ts.mtimeLOC - tref + 1; % to match doy
doy = sample_list(:,6);

% Plot
figure(), set(gcf,'units','centimeters','position',[2 2 50 30])
stairs(tnew,par_ts.data,'-k'), hold on
stairs(doy,PAR2PLOT(:,7),'-r','linewidth',2) % natural day irradiance
stairs(doy - 1 + tsamLOC,PAR2PLOT(:,3),'-b','linewidth',2) % 24h prior to sampling irradiance
xx = [0/48 2/48];
for id = 1:size(sample_list,1)
%     if mod(id,2)
%         tsamLOC = 10.5/24;
%     else
%         tsamLOC = 14.5/24;
%     end
    plot(xx + doy(id)+ noonLOC,PAR2PLOT(id,9)*ones(1,2),'-r','linewidth',5)
    plot(xx + doy(id)+ tsamLOC,PAR2PLOT(id,1)*ones(1,2),'-g','linewidth',5)
    plot(xx + doy(id)+ tsamLOC,PAR2PLOT(id,2)*ones(1,2),'-b','linewidth',5)
    plot(doy(id) - [2 0] + tsamLOC,PAR2PLOT(id,5)*ones(1,2),':o','color',[.3 .3 .8]) % 24h prior to sampling irradiance
end
% plot(doy+real_noonLOC,PAR2PLOT(:,1),'*b')
% plot(doy+tsam_LOC,PAR2PLOT(:,1),'*b')
% plot(doy+tsam_LOC,PAR2PLOT(:,2),'*g')

xlim([tmin_plot tmax_plot])
xlabel('Day of year','fontsize',18)
ylabel('PAR irradiance (µmol m^{-2} s^{-1})','fontsize',18)
title('Ice camp, 1 to 10 May 2016','fontsize',22)
set(gca,'fontsize',14)