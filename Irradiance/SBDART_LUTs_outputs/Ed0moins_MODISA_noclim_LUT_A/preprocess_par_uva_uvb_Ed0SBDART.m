% IMPORT PAR DATA FROM GE ICE CAMP 2016, CALCULATE DAILY PAR AND ALSO PAR
% OF THE 2 DAYS PRIOR TO SAMPLING
clc, clear

%% General settings

% ------------------------------- SELECT -------------------------------

% With LUT A, Ed0moins
datadir = 'Ed0moins_MODISA_noclim_LUT_A';
prodname = 'Edmoins';

% % With LUT B, Ed0moins
% datadir = 'Ed0moins_MODISA_noclim_LUT_B';
% prodname = 'Edmoins';
%
% % With LUT B, Ed0Plus
% datadir = 'Ed0moins_MODISA_noclim_LUT_A';
% prodname = 'EdPlus';

% ------------------------------- SELECT -------------------------------

genpath = '/Users/martigalitapias/Desktop/GreenEdge/Irradiance';
datapath = sprintf('%s/%s',genpath,datadir);
stn_log = 'input.noclim.IceCamp20152016.txt';
sample_list = dlmread(sprintf('%s/%s',genpath,stn_log));
year = 2016;
specres = 5;
wl = 290:specres:700;

%% Create matrix with correlative SBDART output

EdSpec = [];
mtimeUTC = [];
file_list = dir(sprintf('%s/*67.480_-63.790.txt',datapath));

if size(sample_list,1) ~= size(file_list,1)
    error('Some files missing')
else
    for id = 1:size(sample_list,1)
        sl = sample_list(id,:);
        filename = sprintf('Ed0_%4i_%03i_%0.3f_%0.3f.txt',sl(5),sl(6),sl(1),sl(2));
        tmp = dlmread(sprintf('%s/%s',datapath,filename));
        EdSpec = [EdSpec; tmp'];
        mtimeUTC = [mtimeUTC; datenum(sl(5),0,sl(6),(0:3:21)',0,0)];
    end
end

% Calculate local time for Qik (is UTC time minus 4h)
mtimeLOC = mtimeUTC - 4/24;

%% Select time period and integrate spectrum
mtimeUTC_vec = datevec(mtimeUTC);
select = mtimeUTC_vec(:,1)==year;
mtimeUTC = mtimeUTC(select);
mtimeUTC_vec = mtimeUTC_vec(select,:);
EdSpec = EdSpec(select,:);
wlpar = wl>=400 & wl<=700;
par = mean(EdSpec(:,wlpar),2)*300;
wluva = wl>=320 & wl<=400;
uva = mean(EdSpec(:,wluva),2)*80;
wluvb = wl>=290 & wl<=320;
uvb = mean(EdSpec(:,wluvb),2)*20;

%% Plot to check
figure(),
x = datevec(mtimeUTC);
x = yearday(x(:,3),x(:,2),x(:,1),x(:,4),x(:,5),x(:,6));
plot(x,par,'--k'), hold on
plot(x,uva*10,'-b')
plot(x,uvb*400,'-r')
xlabel(sprintf('Day of year (%i)',year))
ylabel('µmol photons m^{-2} s^{-1}')
xlim([110 125])
% xlim([90 230])
legend('PAR','UVA*10','UVB*400');

%% Save data in mat format
par_ts.mtimeLOC = mtimeLOC;
par_ts.mtimeUTC = mtimeUTC;
par_ts.data = par;
par_ts.dataUVA = uva;
par_ts.dataUVB = uvb;
units = 'micromol m-2 s-1';
save(sprintf('%s/%s_PAR_UVA_UVB.mat',genpath,datadir),'par_ts','units')
