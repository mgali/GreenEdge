% MATCH IRRADIANCE AND ICE DATA
clc
clear

%% Load ice data for day of sampling and previous 2 days
icefile = xlsread('~/Desktop/GreenEdge/Ice_Snow/ice_concentration_ice_camp_2015_2016.xlsx','A2:D306');

truedate1 = datenum(2015,3,1,0,0,0);
timelag1 = (icefile(1,2)) - truedate1;
truedate6 = datenum(2015,3,6,0,0,0);
timelag6 = (icefile(6,2)) - truedate6;
if timelag1==timelag6
    icefile(:,2) = icefile(:,2) - timelag1;
else
    error('check dates again!')
end
dv = datevec(icefile(:,2));
doy = yearday(dv(:,3),dv(:,2),dv(:,1),0,0,0);


%% Match by day

stnlist = dlmread('~/Desktop/GreenEdge/Irradiance/input.noclim.IceCamp20152016.txt');
iceout = nan(size(stnlist,1),1);

for j = 1:size(stnlist,1)
    imatch = dv(:,1) == stnlist(j,5) & doy == stnlist(j,6);
    if sum(imatch) == 1
        iceout(j) = icefile(imatch,1)/100;
    elseif sum(imatch) > 1
        error('CHECK DATA AND CODE!')
    else
        iceout(j) = 0.05;
    end
end
stnlist = [stnlist iceout];


%% Plot to check and correct any missing value ad hoc

figure(),plot(stnlist(:,6),stnlist(:,7),'+');

stnlist(7,7) = nanmean(stnlist([6 8],7));

figure(),plot(stnlist(:,6),stnlist(:,7),'+');

%% Write out

dlmwrite('input.noclim.IceCamp20152016.withIce.txt',stnlist,'delimiter',' ','precision','%0.4f');