% IMPORT PAR DATA FROM GE ICE CAMP 2016, CALCULATE DAILY PAR AND ALSO PAR
% OF THE 2 DAYS PRIOR TO SAMPLING
clc, clear

A = xlsread('~/Desktop/GreenEdge/Irradiance/PAR_insitu/GE2016/PAR_GreenEdge_RAD_2016.xlsx','A2:B59041');
K = xlsread('~/Desktop/GreenEdge/Irradiance/PAR_insitu/GE2016/PAR_GreenEdge_RAD_2016.xlsx','D2:D59041');
LW = xlsread('~/Desktop/GreenEdge/Irradiance/PAR_insitu/GE2016/PAR_GreenEdge_RAD_2016.xlsx','F2:F59041');

%% Corrections

% It seems matlab gets the date with -1904 years delay, so let's correct it
dummy = datevec(A(:,1));
dummy(:,1) = dummy(:,1) + 1904; % add 1904 years to get 2016
mtimeUTC = datenum(dummy);

%% Plot to check
figure(),
x = datevec(mtimeUTC);
x = yearday(x(:,3),x(:,2),x(:,1),x(:,4),x(:,5),x(:,6));
plot(x,A(:,2),'-k'), hold on
plot(x,K(:,1),'-r'), 
xlabel('Day of year (2016)'), ylabel('Irradiance')
% xlim([155 200])
xlim([155 160])
legend('PAR, µmol photons m^{-2} s^{-1}','TOT SW W m^{-2}')

%% Save data in mat format

% par_ts.mtimeUTC = mtimeUTC;
% par_ts.data = A(:,2);
% parunits = 'micromol m-2 s-1';
% save('PAR_Data_GE2016-Amundsen_cyber.mat','par_ts','parunits')

%% Calculate conversion factor K to PAR

convfact = 2.77e18; % quanta s-1 W-1 % Morel & Smith 1977, in Kirk 2013 (page 5)
avogadro = 6.022140857e23; % mol-1
PARw = A(:,2)*avogadro./(1e6*convfact);
dPARw = smooth(PARw,60);
dK = smooth(K,60);
fract = dPARw./dK;
quantile(fract,[.1 .25 .5 .75 .9])

lat = 69;
lon = -60;
[azim,elev] = SolarAzEl(mtimeUTC,lat,lon,0);

%%
fract = dPARw./dK;
figure(8),clf
plot(90-elev,fract,'.k'), hold on
fract(fract<quantile(fract,.025) | fract>quantile(fract,.975)) = nan;
figure(8),plot(90-elev,fract,'.b')
ylim([0 1])
xlabel('SZA'), ylabel('PAR/TOT')
