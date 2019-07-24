% INPUT CALCULATION OF PAR DAILY METRICS AT ICE CAMPS. DAY BY DAY INSTEAD
% OF BASED ON SAMPLES LIST

clear
clc

%% Load samples list

year = 2015;

load(sprintf('samples.IceCamp%04i.mat',year))
if year == 2015
    samples = samplesIC2015;
    header = headerIC2015;
elseif year == 2016
    samples = samplesIC2016;
    header = headerIC2016;
    note = noteIC2016;
end

% Doy vector
doy0 = yearday(3,4,year,0,0,0);
doyf = yearday(31,7,year,0,0,0);
doy = doy0:doyf;

% Mean sampling time
meanUTC = nanmean(samples(:,5));
medianUTC = nanmedian(samples(:,5));

% Coords
lat = mode(samples(:,6));
lon = mode(samples(:,7));

% Preallocate
samplesIC = nan(length(doy),length(header));
samplesIC(:,2) = doy;

% Fill
for j = 1:length(doy)
    
    tmp = samples(samples(:,2) == doy(j),:);
    
    if isempty(tmp)
        samplesIC(j,1:7) = [year doy(j) NaN NaN medianUTC lat lon];
    else
        samplesIC(j,:) = nanmean(tmp);
    end
end


%% Save txt
dlmwrite(sprintf('samplesAllDays.IceCamp%04i.txt',year),samplesIC,'delimiter',' ','precision','%0.4f');

%% Save mat

if year == 2015
    headerIC2015 = header;
    samplesIC2015 = samplesIC;
    samplesIC2015(samplesIC2015==-999) = nan;
    save(sprintf('samplesAllDays.IceCamp%04i.mat',year),'samplesIC2015','headerIC2015')
elseif year == 2016
    headerIC2016 = header;
    noteIC2016 = note;
    samplesIC2016 = samplesIC;
    samplesIC2016(samplesIC2016==-999) = nan;
    save(sprintf('samplesAllDays.IceCamp%04i.mat',year),'samplesIC2016','headerIC2016','noteIC2016')
end

