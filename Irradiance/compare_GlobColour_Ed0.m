% Test mean_par_period.m on 2016 ice camp PAR data, plot
clc, clear

%% General definitions

tsearch_h = 3.1;
pXh = 48;
UTC_diff = 4;
lon = -63.78953;

specres = 5;
wl = 290:specres:700;
wlpar = wl>=400 & wl<=700;

% -------------------------- SELECT IN SITU DATASET -----------------------

icecamp = 0;
year = 2016;

if icecamp
    % Ice camp 2016
    load ~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2016-ICECAMP_cyber.mat
    stn_list = dlmread('input.noclim.IceCamp20152016.txt');
    load samples.IceCamp2016.mat
    stn_list(stn_list(:,5)~=year,:) = [];
    samples = samplesIC2016;
    samples(:,5) = 24*samples(:,5);
    header_samples = headerIC2016;
    tm_insitu = 1;
    path_GlobColour = '~/Desktop/GreenEdge/Irradiance/matchups_GlobColour/Matchup_Globcolour_IceCamp20152016_4Matlab.txt';
else
    % Amundsen 2016
    load ~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2016-Amundsen_cyber.mat
    stn_list = dlmread('input.noclim.GE2016.txt');
    load samples.GE2016.castsOnly.mat
    stn_list(stn_list(:,5)~=year,:) = [];
    samples = samplesGE2016castsOnly;
    samples(:,5) = 24*samples(:,5);
    header_samples = headerGE2016;
    tm_insitu = 1/60;
    path_GlobColour = '~/Desktop/GreenEdge/Irradiance/matchups_GlobColour/Matchup_Globcolour_noclim.GE2016_4Matlab.txt';
end

% ---------------- Modify samples list to include all days ----------------
STMP = samples;
samples = nan(size(stn_list,1),7);
samples(:,[1 2 6 7]) = stn_list(:,[5 6 1 2]);
samples(:,5) = nanmean(STMP(:,5))*ones(size(samples,1),1)/24;
samples(samples(:,2)<94,:) = 94;
samples(samples(:,2)>227,:) = 227;

% --------------------------- SELECT SBDART DATASET -----------------------

% % With LUT A, Ed0moins
% datadir = 'Ed0moins_MODISA_noclim_LUT_A';
% prodname = 'Edmoins';

% % With LUT B, Ed0moins
% datadir = 'Ed0moins_MODISA_noclim_LUT_B';
% prodname = 'Edmoins';

% With LUT B, Ed0Plus
datadir = 'Ed0Plus_MODISA_noclim_LUT_B';
prodname = 'EdPlus';

% ---------------------------------- END SELECT ---------------------------

genpath = '/Users/martigalitapias/Desktop/GreenEdge/Irradiance';
datapath = sprintf('%s/%s',genpath,datadir);


%% Load GlobColour PAR data

[head,DATA] = f_read_txt_dummy(path_GlobColour,1,1,' ',1,1);
h_GlobColour = head';
PAR0_GlobColour = nan(size(DATA,1),length(head));
for k = [1 3:10]
    for j = 1:size(DATA,1)
        PAR0_GlobColour(j,k) = str2double(DATA{j,k});
    end
end
for j = 1:size(DATA,1)
    PAR0_GlobColour(j,2) = datenum(DATA{j,2});
end

% Append datevec format
DV = datevec(PAR0_GlobColour(:,2));
PAR0_GlobColour(:,11:13) = DV(:,1:3);
PAR0_GlobColour(:,14) = yearday(DV(:,3),DV(:,2),DV(:,1),DV(:,4),DV(:,5),DV(:,6));
h_GlobColour = [h_GlobColour; ({'year' 'month' 'day' 'doy'})'];


%% Preallocate and fill in situ and SBDART PAR data

PAR0_insitu = nan(size(stn_list,1),21);
PAR0_sbdart = nan(size(stn_list,1),21);
tsamplUTC = nan(size(stn_list,1),1);
real_noonLOC = nan(size(stn_list,1),1);

for id = 1:size(stn_list,1)
    
    sl = stn_list(id,:);
    
    if icecamp
        tmp = samples(samples(:,1)==year & samples(:,2)==sl(6),:);
    else
        if size(samples,1) == size(stn_list,1), tmp = samples(id,:);
        else error('Number of casts in samples and stn_list must match in the Amundsen data set')
        end
    end
    UTC_sample = nanmean(tmp(:,5));
    if isnan(UTC_sample), UTC_sample = nanmean(samples(:,5)); end
    
    if ~isempty(tmp) && ~isnan(UTC_sample)
        
        tsamplUTC(id) = datenum(sl(5),sl(3),sl(4),UTC_sample,0,0);
        real_diff = -sl(2)/15;
        real_noonLOC(id) = 0.5 + (real_diff-UTC_diff)/24;
        
        % In situ PAR data
        PARTMP = mean_par_period(par_ts, tsamplUTC(id), tsearch_h, pXh, nanmean(tmp(:,7)), UTC_diff, tm_insitu);
        PAR0_insitu(id,:) = (cell2mat(struct2cell(PARTMP)))';
        
        % Satellite SBDART-derived data
        backsearch = [-2 -1 0 1];
        EdSpec = [];
        mtimeUTC = [];
        for k = 1:length(backsearch)
            tmpEdSpec = dlmread(sprintf('%s/Ed0_%4i_%03i_%0.3f_%0.3f.txt',datapath,sl(5),sl(6)+backsearch(k),sl(1),sl(2)));
            EdSpec = [EdSpec; tmpEdSpec'];
            mtimeUTC = [mtimeUTC; datenum(sl(5),0,sl(6)+backsearch(k),(0:3:21)',0,0)];
        end
        par_ts_sbdart.data = mean(EdSpec(:,wlpar),2)*300;
        par_ts_sbdart.mtimeUTC = mtimeUTC;
        par_ts_sbdart.mtimeLOC = mtimeUTC - UTC_diff/24;
        PARTMP = mean_par_period(par_ts_sbdart, tsamplUTC(id), tsearch_h, pXh, nanmean(tmp(:,7)), UTC_diff, 3);
        PAR0_sbdart(id,:) = (cell2mat(struct2cell(PARTMP)))';
    end
    
    if icecamp
        % Find matching index for GlobColour matchups
        [dist_id, az] = distance(PAR0_GlobColour(:,3),PAR0_GlobColour(:,4),sl(1),sl(2));
        matchlog = PAR0_GlobColour(:,11)==year & PAR0_GlobColour(:,14)==sl(6) & dist_id*60*1.852 < 10;
        if sum(matchlog)
            indmatch{id} = find(matchlog == 1);
        end
    end
    
end

tsamLOC = tsamplUTC - UTC_diff;
real_noonUTC = real_noonLOC + UTC_diff;

% Reorder PAR_GlobColour
if icecamp
    GTMP = PAR0_GlobColour;
    PAR0_GlobColour = nan(size(PAR0_GlobColour));
    for j = 1:size(stn_list,1)
        ij = indmatch{j};
        if length(ij)>1, ij, end
        PAR0_GlobColour(j,:) = GTMP(min(ij),:);
    end
end
if size(PAR0_GlobColour,1) > size(stn_list,1)
    PAR0_GlobColour(size(stn_list,1)+1:end,:) = [];
end

%% --------- Plot original time series along with extracted values ---------

% Assign variables
h_mean_par_period_output = ({'tsam' 'tsami' 'tsam3h' 'Ntsam3h' 'tsam1h'...
    'Ntsam1h' 'p24h' 'Np24h' 'p48h' 'Np48h'...
    'p24hmax' 'day' 'Nday' 'noon' 'noon1h'...
    'Nnoon1h' 'noon3h' 'Nnoon3h' 'daymax' 'dayUTC' 'NdayUTC'})';


% vnames = ({'tsam' 'tsami' 'tsam3h'})'; vcolors = [0 0 1; 0 1 0; 1 0 0];
vnames = ({'p24h' 'p48h' 'day'})'; vcolors = [0 0 1; 0 1 0; 1 0 0];
% vnames = ({'tsam' 'tsam3h'})'; vcolors = [0 0 1; 1 0 0];
sc = 'log'; % 'lin' or 'log'
X = nan(size(stn_list,1),length(vnames));
Y = nan(size(X));

%% Comparison with multiple scatterplots

X = nan(size(stn_list,1),2);
X(:,1) = PAR0_insitu(:,strcmp('dayUTC',h_mean_par_period_output));
X(:,2) = PAR0_sbdart(:,strcmp('dayUTC',h_mean_par_period_output));
xl = {'In situ' 'SBDART-Aqua' 'In situ'};

yconvfact = 1e6/(3600*24); % mol quanta m-2 d-1 to µmol quanta m-2 s-1
Y = PAR0_GlobColour(:,5:10)*yconvfact;

if ~icecamp
    load ~/Desktop/GreenEdge/Irradiance/CASTS_ICE_matched_24-Mar-2017.mat
    iceday = nanmean(OUT(:,(end-1):end),2);
end
yl = h_GlobColour(5:10);

xymin = 0;
xymax = nanmax([X(:); Y(:)]);
icecut = 0.75;

figure(321), clf
cc = 0;
for k = 1:size(X,2)
    for j = 1:size(Y,2)
        cc = cc + 1;
        % subplot(size(Y,2),size(X,2),cc) % jXk plot
        subplot(size(X,2),size(Y,2),cc) % kXj plot
        scatter(X(:,k),Y(:,j),30,'b','filled'), hold on
        plot([xymin xymax],[xymin xymax],'-r')
        plot([xymin xymax],1.5*([xymin xymax]),':r')
        plot([xymin xymax],0.5*([xymin xymax]),':r')
        axis([xymin xymax xymin xymax])
        axis square
        xlabel(xl{k},'interpreter','none','fontsize',12)
        ylabel(yl{j},'interpreter','none','fontsize',12)
        set(gca,'fontsize',12,'xtick',0:200:600,'ytick',0:200:600)
        grid on
        
        s = f_skill_stats(X(:,k),Y(:,j),'lin','stats4diagramsOFF');
        if ~isempty(s)
            text(20,700,sprintf('r = %0.2f',s.r),'fontsize',10,'color','r')
            text(20,640,sprintf('RMSE = %0.0f',s.rms),'fontsize',10,'color','r')
            text(20,570,sprintf('Bias = %0.0f%%',100*s.bias/nanmean(X(:,k))),'fontsize',10,'color','r')
            % text(20,510,sprintf('Slope = %0.2f',s.slope),'fontsize',10,'color','r')
        end
        
        if ~icecamp
            scatter(X(iceday>icecut,k),Y(iceday>icecut,j),30,'c','filled')
        end
        if cc == 1
            scatter(1000,1000,30,'c','filled')
            legend('All','1:1','1.5:1','0.5:1',sprintf('Ice>%0.2f',icecut),'location','southwest') % northwest
        end
        
    end
end


%% Ratio, for icecamp

% if icecamp
%     figure(112), clf, hold on
%     set(gcf,'units','centimeters','position',[3 25 30 20])
%     xax = stn_list(:,6);
%     bar(xax,Y./X)
%     xlim([nanmin(xax)-2 nanmax(xax)+2])
%     legend(vnames,'location','northwest')
%     plot([nanmin(xax)-2 nanmax(xax)+2],[1 1],'-k')
%     xlabel(sprintf('Day of year %4i',year),'fontsize',20)
%     ylabel('Ratio SBDART-Aqua / In situ','fontsize',20)
%     set(gca,'fontsize',16)
% end
