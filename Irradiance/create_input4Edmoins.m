% INPUT FOR RADIATIVE TRANSFER CODE (SBDART) AS USED FOR GREEN EDGE CRUISE
% MODIS-A DATA USED (SINCE ISCCP IS NOT AVAILABLE AFTER 2009)
% CLIMATOLOGY MODE HAS NOT BEEN IMPLEMENTED FOR MODIS-A ATMOSPHERE DATA

clear
clc

modays = [31 28 31 30 31 30 31 31 30 31 30 31];
modaysleap = [31 29 31 30 31 30 31 31 30 31 30 31];
mocum = cumsum(modays); mocum = [0 mocum(1:(end-1))];
mocumleap = cumsum(modaysleap); mocumleap = [0 mocumleap(1:(end-1))];

%% Define periods and coordinates, create stations-dates list

% Years and months
years = [2015 2016];
months = 4:8;

% Create years and months
mdy = [];
for iy = years
    for im = months
        if ~mod(iy,4), md = modaysleap(im); else md = modays(im); end
        if im==8, md = 15; end % Trick to run code only for first 15 days of August
        for id = 1:md
            mdy = [mdy; [im id iy]];
        end
    end
end
month = mdy(:,1);
day = mdy(:,2);
year = mdy(:,3);

%% Calculate doy
doy = nan(size(mdy,1),1);
isleap = ~mod(year,4);
doy(isleap) = (mocumleap(month(isleap)))' + day(isleap);
isnotleap = ~~mod(year,4);
doy(isnotleap) = (mocum(month(isnotleap)))' + day(isnotleap);

% % Same thing done in loop. Just to check that logical indexing is well done
% doyloop = nan(size(day));
% for i = 1:length(doyloop)
%     if ~mod(year(i),4)
%         doyloop(i) = mocumleap(month(i)) + day(i);
%     else
%         doyloop(i) = mocum(month(i)) + day(i);
%     end
% end

%% Add sampling coordinates

latIC = 67 + 28.784/60; % ice camp lat
lonIC = -(63 + 47.372/60);  % ice camp lon

lat = latIC*ones(size(mdy,1),1);
lon = lonIC*ones(size(mdy,1),1);

%% Create stations-dates list. Write text file

OUT = [lat lon month day year doy];
% filename = 'input.noclim.IceCamp20152016';
% dlmwrite(strcat(filename,'.txt'),OUT,'delimiter',' ','precision','%0.3f');


%% ------------------------------------------------------------------------
% Load ship data

% ! head -1 station_data/LogBookScience_GreenEdge_leg1a.csv.txt
shipheader = {'day' 'month' 'year' 'H' 'M' 'dayUTC' 'monthUTC' 'yearUTC' 'HUTC' 'MUTC' 'latNdegr'...
    'latNmin' 'lonWdegr' 'lonWmin' 'cap' 'cast' 'zbot' 'Wdir' 'Wspeed' 'Tair' 'Twater' 'Patm' 'RH' 'Ice'};
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

GOUT = [lat lon month day year doy];
filename = 'input.noclim.GE2016';
dlmwrite(strcat(filename,'.txt'),GOUT,'delimiter',' ','precision','%0.3f');

%% Append IceCamp and GE cruise
ALLOUT = [OUT; GOUT];
filename = 'input.noclim.GEALL';
dlmwrite(strcat(filename,'.txt'),ALLOUT,'delimiter',' ','precision','%0.3f');
