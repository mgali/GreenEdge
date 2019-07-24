% INPUT FOR RADIATIVE TRANSFER CODE (SBDART) AS USED FOR GREEN EDGE CRUISE
% MODIS-A DATA USED (SINCE ISCCP IS NOT AVAILABLE AFTER 2009)
% CLIMATOLOGY MODE HAS NOT BEEN IMPLEMENTED FOR MODIS-A ATMOSPHERE DATA

clear
clc

%% General definitions
modays = [31 28 31 30 31 30 31 31 30 31 30 31];
modaysleap = [31 29 31 30 31 30 31 31 30 31 30 31];
mocum = cumsum(modays); mocum = [0 mocum(1:(end-1))];
mocumleap = cumsum(modaysleap); mocumleap = [0 mocumleap(1:(end-1))];
UTCdiff = 4/24;
latIC = 67 + 28.784/60;
lonIC = -63.78953;

%% General settings
dir_stndata = '~/Desktop/GreenEdge/Irradiance/station_data';

%% IceCamp 2015 (no PAR sensor data available)

fname2015 = 'LogElectroniqueOperationsGreenEdge2015_1_Marti.xlsx';
dummyIC2015 = xlsread(sprintf('%s/%s',dir_stndata,fname2015),'B2:N694');
samplesIC2015 = [2015*ones(size(dummyIC2015,1),1) dummyIC2015(:,[1 13 12 9]) latIC*ones(size(dummyIC2015,1),1) lonIC*ones(size(dummyIC2015,1),1)];
samplesIC2015(:,5) = samplesIC2015(:,5) + UTCdiff;

% Save txt
dlmwrite('samples.IceCamp2015.txt',samplesIC2015,'delimiter',' ','precision','%0.4f');

% Save mat
headerIC2015 = ({'year' 'doy' 'bottle' 'depth' 'UTC_decimal' 'lat' 'lon'})';
save('samples.IceCamp2015.mat','samplesIC2015','headerIC2015');

%% IceCamp 2016 (hourly PAR sensor data available during certain period)

fname2016 = 'GE21016-Water_sample_logbook.xlsx';
dummyIC2016 = xlsread(sprintf('%s/%s',dir_stndata,fname2016),'A2:E229');
samplesIC2016 = [2016*ones(size(dummyIC2016,1),1) dummyIC2016(:,2:5) latIC*ones(size(dummyIC2016,1),1) lonIC*ones(size(dummyIC2016,1),1)];
samplesIC2016(:,5) = samplesIC2016(:,5) + UTCdiff;

% Save txt
dlmwrite('samples.IceCamp2016.txt',samplesIC2016,'delimiter',' ','precision','%0.4f');

% Save mat
headerIC2016 = ({'year' 'doy' 'bottle' 'depth' 'UTC_decimal' 'lat' 'lon'})';
noteIC2016 = {'Depth=0 means under-ice sample'};
save('samples.IceCamp2016.mat','samplesIC2016','headerIC2016','noteIC2016');


%% ------------------------------------------------------------------------
% Load ship data

% ! head -1 station_data/LogBookScience_GreenEdge_leg1a.csv.txt
shipheader = ({'day' 'month' 'year' 'H' 'M' 'dayUTC' 'monthUTC' 'yearUTC' 'HUTC' 'MUTC' 'latNdegr'...
    'latNmin' 'lonWdegr' 'lonWmin' 'cap' 'cast' 'zbot' 'Wdir' 'Wspeed' 'Tair' 'Twater' 'Patm' 'RH' 'Ice'})';
leg1a = dlmread('~/Desktop/GreenEdge/Irradiance/station_data/LogBookScience_GreenEdge_leg1a.csv.nohead.txt');
leg1b = dlmread('~/Desktop/GreenEdge/Irradiance/station_data/LogBookScience_GreenEdge_leg1b.csv.nohead.txt');
leg1b(:,end) = [];
G = [leg1a; leg1b];

lat = G(:,11) + G(:,12)/60;
lon = -(G(:,13) + G(:,14)/60);
month = G(:,7);
year = G(:,6);
day = G(:,8);

doy = nan(size(year,1),1);
isleap = ~mod(year,4);
doy(isleap) = (mocumleap(month(isleap)))' + day(isleap);
isnotleap = ~~mod(year,4);
doy(isnotleap) = (mocum(month(isnotleap)))' + day(isnotleap);

UTC_decimal = (G(:,9) + G(:,10)/60)/24;

samplesGE2016castsOnly = [year doy -999*ones(size(year,1),2) UTC_decimal lat lon];

% Save txt
dlmwrite('samples.GE2016.castsOnly.txt',samplesGE2016castsOnly,'delimiter',' ','precision','%0.4f');

% Save mat
headerGE2016 = ({'year' 'doy' 'bottle' 'depth' 'UTC_decimal' 'lat' 'lon'})';
noteGE2016 = {'Cast data only, depths and bottle number not included so far'};
samplesGE2016castsOnly(samplesGE2016castsOnly==-999) = nan;
save('samples.GE2016.castsOnly.mat','samplesGE2016castsOnly','headerGE2016','noteGE2016')

