% IMPORT PAR DATA FROM GE ICE CAMP 2015, CALCULATE DAILY PAR AND ALSO PAR
% OF THE 2 DAYS PRIOR TO SAMPLING
clc, clear

datapath = '~/Desktop/GreenEdge/Irradiance/PAR_insitu/IC2015';
filename = 'RadiationData-GreenEdge-1';

dtpre = xlsread(sprintf('%s/%s.xlsx',datapath,filename),'A2:A51143'); % date and time
edpre = xlsread(sprintf('%s/%s.xlsx',datapath,filename),'E2:E51143'); % downwelling shortwave irradiance
apre = xlsread(sprintf('%s/%s.xlsx',datapath,filename),'M2:M51143'); % albedo

%% It seems matlab gets the date with a one-day and 1900 years delay, so
% let's correct it
dummy = datevec(dtpre);
dummy(:,1) = dummy(:,1) + 1900; % add 1900 years to get 2016
dummy2 = datenum(dummy);
mtimeUTC = dummy2 - 1;

%% Correct baseline of irradiance
edpre(edpre<0) = 0;

%% Plot to check irradiance
figure(),
x = datevec(mtimeUTC);
doy = yearday(x(:,3),x(:,2),x(:,1),x(:,4),x(:,5),x(:,6));
plot(doy,edpre,'.r')
xlabel('Day of year (2016)'), ylabel('Ed shortwave, W m^{-2}')

%% Convert to quanta
convfact = 2.77e18; % quanta s-1 W-1 % Morel & Smith 1977, in Kirk 2013 (page 5)
avogadro = 6.022140857e23; % mol-1

% Min and Max PAR/TOTAL_SW fraction
parfract_min = 0.42; % at SZA > 80;
parfract_max = 0.45; % at SZA < 50;

% Load coordinates
load('~/Desktop/GreenEdge/Irradiance/samples.IceCamp2015.mat')
lat = samplesIC2015(1,6);
lon = samplesIC2015(1,7);
[azim,elev] = SolarAzEl(mtimeUTC,lat,lon,0);
sza = 90 - elev;

% % Estimate PAR fraction as function of SZA
% parfract = nan(size(doy));
% parfract(sza>80) = parfract_min;
% parfract(sza<50) = parfract_max;
% subset = sza<=80 & sza >= 50;
% parfract(subset) = interp1([50; 80],[parfract_max; parfract_min],sza(subset));
parfract = 0.43;

%% Plot
figure(2),
plot(doy,sza,'.r'), hold on
plot(doy,elev,'.b')
plot(doy,100*parfract,'.k')
xlabel('Day of year (2016) ','fontsize',20), ylabel('Angle (degrees) ','fontsize',20)
legend('Solar zenith angle ','Solar elevation ','PAR/PYR fraction ')
xlim([145 150])
plot([145 150],[10 10],'-','color',[.7 .7 .7])
plot([145 150],[40 40],'-','color',[.7 .7 .7])
plot([145 150],[50 50],'-','color',[.7 .7 .7])
plot([145 150],[80 80],'-','color',[.7 .7 .7])
set(gca,'fontsize',16)


%% Correct
% parfract = 0.43;
par = edpre.*parfract*convfact*1e6/avogadro;

%% Plot to check irradiance
figure(),
plot(doy,par,'.b')
xlabel('Day of year (2016)'), ylabel('PAR, µmol photons m^{-2} s^{-1}')

%% Save data in mat format
par_ts.mtimeUTC = mtimeUTC;
par_ts.data = par;
parunits = 'micromol m-2 s-1';
save(sprintf('%s/PAR_Data_GE2015-ICECAMP_cyber.mat',datapath),'par_ts','parunits')
