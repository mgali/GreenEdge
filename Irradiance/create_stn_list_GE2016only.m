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
stn = ([[nan nan 100 nan 101:102 nan 104 103 105:107 nan nan 108 109 110 nan nan 111:115 nan nan],... % casts 1-26, n=26
    [200 nan*ones(1,3) 202:203 nan*ones(1,3) 205:207 207 nan nan 208:211 nan],... % casts 27:46, n=20
    [300 nan*ones(1,3) 301:302 nan 304:306 nan*ones(1,3) 307:308 nan 311 313 312 nan nan 314 nan 316:318 nan nan 319:320 nan 322:323 nan*ones(1,3) 325 nan nan],... % casts 47-85, n=39
    [400:403 403 nan 404:409 409 nan 410:415 nan*ones(1,4) 418 nan nan 416:417 419:420],... % casts 86-116, n=31
    [500:505 507 507 nan 506.5 506 508:509 nan 511:513 512 nan nan 513:519 nan nan 521],... % casts 117-146, n=30
    [600 600 nan 601:605 605 605 nan 606 615 615 nan 614 613 612 611 610 604.5 604.5 607:609 616 617 nan],... % casts 147-174, n=28
    [703 nan nan 700 701:702 704 705 707 nan*ones(1,3) 709 710 nan*ones(1,3) 713 713 713 712 714:719 nan nan]])'; % casts 175-203, n=29

doy = nan(size(year,1),1);
isleap = ~mod(year,4);
doy(isleap) = (mocumleap(month(isleap)))' + day(isleap);
isnotleap = ~~mod(year,4);
doy(isnotleap) = (mocum(month(isnotleap)))' + day(isnotleap);

UTC_decimal = (G(:,9) + G(:,10)/60)/24;

samplesGE2016castsOnly = [year doy -999*ones(size(year,1),2) UTC_decimal lat lon stn];

% Save txt
dlmwrite('samples.GE2016.castsOnly.txt',samplesGE2016castsOnly,'delimiter',' ','precision','%0.4f');

% Save mat
headerGE2016 = ({'year' 'doy' 'bottle' 'depth' 'UTC_decimal' 'lat' 'lon' 'stn'})';
noteGE2016 = {'Cast data only, depths and bottle number not included so far'};
samplesGE2016castsOnly(samplesGE2016castsOnly==-999) = nan;
save('samples.GE2016only.castsOnly.mat','samplesGE2016castsOnly','headerGE2016','noteGE2016')

