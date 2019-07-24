% IMPORT PAR DATA FROM GE ICE CAMP 2016, CALCULATE DAILY PAR AND ALSO PAR
% OF THE 2 DAYS PRIOR TO SAMPLING
clc, clear

datapath = '~/Desktop/GreenEdge/Irradiance/PAR_insitu/IC2015';
filename = 'GE15_OCR_Ed.dat';
fid = fopen(sprintf('%s/%s',datapath,filename));
datacell = textscan(fid,'%f%f%s%f%f%f%f%f%f%f%f%f%f','delimiter',',','HeaderLines',1);
fclose(fid);

par_ts.lat = datacell{1};
par_ts.lon = datacell{2};
par_ts.mtimeUTC = nan(size(par_ts.lat));

%% Manage date and time
dtpre = datacell{3};
dummy = regexp(dtpre,' ','split');

for j = 1:length(dummy)
    tmp = dummy{j,1};
    tmp1 = regexp(tmp{1,1},'/','split');
    tmp2 = regexp(tmp{1,2},':','split');
    par_ts.mtimeUTC(j) = datenum(str2double(tmp1{3}),str2double(tmp1{2}),str2double(tmp1{1}),str2double(tmp2{1}),str2double(tmp2{2}),str2double(tmp2{3}));
end
    
%% Spectral irradiance
wl = [406.42 434.48 442.33 489.70 508.43 560.22 624.20];
EDspec = [datacell{4} datacell{5} datacell{6} datacell{7} datacell{8} datacell{9} datacell{10}];


%% Plot to check
figure(),
x = datevec(par_ts.mtimeUTC);
x = yearday(x(:,3),x(:,2),x(:,1),x(:,4),x(:,5),x(:,6));
plot(x,EDspec,'.','markersize',4)
xlabel('Day of year (2016)'), ylabel('Ed(wl) µW cm^{-2} s^{-1}')
legend(num2str(wl'))

% Don't save, cannot be used for PAR product because the record is not
% continuous