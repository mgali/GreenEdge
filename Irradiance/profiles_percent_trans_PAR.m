% PROCESS CLEAN AND LOESS-SMOOTHED COPS PAR PROFILES (GUISLAIN OUTPUT)
% AND EXTRACT 10 CM PERCENT IRRADIANCE WITH RESPECT TO SURFACE REFERENCE
clear, clc

%% Define dataset and paths

genpath = '~/Desktop/GreenEdge/Irradiance/COPS_profiles';
datasets = {'GE2015.Ice.Camp';
    'GE2016.Ice.Camp';
    'GE2016.Amundsen.Ice.Stations';
    'GE2016.Amundsen.Open.Water'};
ds = 1;
if ds == 1, year = 2015; else year = 2016; end
dpath = sprintf('%s/%s',genpath,datasets{1});

% % Example folder/file names: quite variable! need flexible code
% 2015163_ICEPRO_150612_1419_C_data_004.tsv.asc
% 20160626_ICEPRO_CAST_009_160626_132004_URC.tsv.asc
% 20160620_COPS_CAST_008_160620_032306_URC.tsv.asc

dirlist = dir(dpath);
dirname = dirlist(4).name;
dummy = regexp(dirname,'_','split');
date = regexp(dummy(1),num2str(year),'split');
date = date{1}{2};
date = char(date);
if length(date) == 3
    doy = str2double(date);
    tm = datenum(year,0,doy,0,0,0);
    tmp = datevec(tm);
    mm = tmp(2);
    dd = tmp(3);
elseif length(date) == 4
    mm = str2double(date(1:2));
    dd = str2double(date(1:2));
    tm = datenum(year,mm,dd,0,0,0);
    doy = yearday(dd,mm,year,0,0,0);
end

%% Load files
tic
% Load PAR
path_PAR = sprintf('%s/%s/PAR.txt',dpath,dirname);
[head,DATA] = f_read_txt_dummy(path_PAR,1,1,'\t',1,1);

% Assign PAR variables
kdate = strcmp(head,'date');
kdepth = strcmp(head,'Depth');
kvalid = strcmp(head,'valid');
kPARd = strcmp(head,'PAR.d');
kPARu = strcmp(head,'PAR.u');
kPAR0 = strcmp(head,'PAR.0');

% Fill PAR
for j = 1:size(DATA,1)
    PAR.date(j,1) = datenum(DATA(j,kdate));
    PAR.z(j,1) = str2double(DATA(j,kdepth));
    PAR.valid(j,1) = str2double(DATA(j,kvalid));
    PAR.d(j,1) = str2double(DATA(j,kPARd));
    PAR.u(j,1) = str2double(DATA(j,kPARu));
    PAR.o(j,1) = str2double(DATA(j,kPAR0));
end

% Load PARfit
path_PARfit = sprintf('%s/%s/PAR.fitted.txt',dpath,dirname);
[head,DATA] = f_read_txt_dummy(path_PARfit,1,1,'\t',1,1);

% Assign PARfit variables
kdepth = strcmp(head,'depth');
kPARd = strcmp(head,'PAR.d.fitted');
kPARu = strcmp(head,'PAR.u.fitted');
kPAR0 = strcmp(head,'PAR.0.fitted');

% Fill PARfit
for j = 1:size(DATA,1)
    PARfit.z(j,1) = str2double(DATA(j,kdepth));
    PARfit.d(j,1) = str2double(DATA(j,kPARd));
    PARfit.u(j,1) = str2double(DATA(j,kPARu));
    PARfit.o(j,1) = str2double(DATA(j,kPAR0));
end

save(sprintf('%s/%s.mat',dpath,dirname),'PAR','PARfit')

toc

%% Plot to examine

figure(31), clf
plot(PAR.d,-PAR.z,'.k'), hold on
plot(PARfit.d,-PARfit.z,'-r')
plot(PAR.u,-PAR.z,'.k')
plot(PARfit.u,-PARfit.z,'-r')
ylim([nanmin(-PAR.z) nanmax(-PAR.z)])
set(gca,'xscale','log')


