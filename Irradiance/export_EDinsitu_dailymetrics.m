% Test mean_par_period.m on 2016 ice camp PAR data, plot
clc, clear


%% ---------------------------------- EDIT --------------------------------

% In situ dataset: 0 is ice camp, 1 is Amundsen
% icecamp = 1; year = 2015;
% icecamp = 1; year = 2016;
icecamp = 0; year = 2016;
% mean_par_period.m settings
tsearch_h = 3.1;
pXh = 48;
UTC_diff = 4;
% SBDART settings
hperiod = 1;
specres = 5;
wl = 290:specres:700;
% =================== Select different spectral ranges =================== 
specrange = 'UVR'; % either PAR, UVR, UVA or UVB
% ========================================================================
if strcmp(specrange,'PAR'), wlmin = 400; wlmax = 700; % PAR
elseif strcmp(specrange,'UVR'), wlmin = 290; wlmax = 400; % UVR
elseif strcmp(specrange,'UVA'), wlmin = 320; wlmax = 400; % UVA
elseif strcmp(specrange,'UVB'), wlmin = 290; wlmax = 320; % UVB
end
wlint = wl>=wlmin & wl<=wlmax;
wlrange = wlmax-wlmin;
hours = (0:hperiod:23)';
% SBDART data version
%   0 is LUT with fixed albedo
%   1 is LUT with variable albedo and simple ice type albedo param
sbdart_data_version = 1;

% Correct bias?
correct_bias = 0;

%% General definitions
if icecamp
    stn_list = dlmread('input.noclim.IceCamp20152016.txt');
    stn_list(stn_list(:,5)~=year,:) = [];
    % load(sprintf('samples.IceCamp%i.mat',year))
    load(sprintf('samplesAllDays.IceCamp%i.mat',year))
    if year == 2015
        %         load ~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2015-ICECAMP_cyber_fract42.mat
        %         load ~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2015-ICECAMP_cyber_fract43.mat
        load ~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2015-ICECAMP_cyber_fractSZA41_44.mat
        samples = samplesIC2015;
        header_samples = headerIC2015;
        tm_insitu = 1/60;
    elseif year == 2016
        load ~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2016-ICECAMP_cyber.mat
        samples = samplesIC2016;
        header_samples = headerIC2016;
        tm_insitu = 1;
    end
    samples(:,5) = 24*samples(:,5);
    backsearch = [-2 -1 0 1];
else
    load ~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2016-Amundsen_cyber.mat
    stn_list = dlmread('input.noclim.GE2016.txt');
    load samples.GE2016.castsOnly.mat
    stn_list(stn_list(:,5)~=year,:) = [];
    samples = samplesGE2016castsOnly;
    samples(:,5) = 24*samples(:,5);
    header_samples = headerGE2016;
    tm_insitu = 1/60;
    backsearch = [-2 -1 0 1];
end

if sbdart_data_version == 0
    % datadir = 'Ed0moins_MODISA_noclim_LUT_A'; % LUT A Ed0moins (Takuvik chain)
    % datadir = 'Ed0moins_MODISA_noclim_LUT_B'; % LUT B Ed0moins (Srikanth +2ozone levels)
    datadir = 'Ed0Plus_MODISA_noclim_LUT_B'; % LUT B Ed0Plus (Srikanth +2ozone levels)
    tm_sbdart = 3; % SBDART temporal resolution in h
elseif sbdart_data_version == 1
    % LUT Ed0plus, variable albedo (SurfAlb folders) and different t resolution
    tm_sbdart = 1; % SBDART temporal resolution in h
    if icecamp
        %     datadir = sprintf('Ed0plus_MODISA_LUT_SurfAlb_v1_%ih_realIce',hperiod);
        %     datadir = sprintf('Ed0plus_MODISA_LUT_SurfAlb_v1_%ih_fakeIce',hperiod);
        datadir = sprintf('Ed0plus_MODISA_LUT_SurfAlb_v1_%ih_consensus',hperiod);
    else
        datadir = sprintf('Ed0plus_MODISA_LUT_SurfAlb_v3_%ih',hperiod);
    end
end

genpath = '/Users/martigalitapias/Desktop/GreenEdge/Irradiance/SBDART_LUTs_outputs';
datapath = sprintf('%s/%s',genpath,datadir);
outpath = '/Users/martigalitapias/Desktop/GreenEdge/Irradiance/PAR0_DAILY_METRICS';

%% Preallocate and fill in situ and SBDART PAR data

% Assign variables
h_mean_par_period_pre = (...
    {'PARsam' 'PARsam_interp' 'PARp1h' 'Np1h' 'PARp3h' 'Np3h' 'PAR3hSAT' 'N3hSAT'...
    'PARp24h' 'Np24h' 'PARp48h' 'Np48h' 'PARp24hmax' 'PARdayLOC' 'NdayLOC' 'PARnoonLOC'...
    'PARnoon1hLOC' 'Nnoon1hLOC' 'PARnoon3hSAT' 'Nnoon3hSAT' 'PARdayLOCmax' 'PARdayUTC' 'NdayUTC'})';

PAR0_insitu = nan(size(stn_list,1),length(h_mean_par_period_pre));
PAR0_sbdart = nan(size(stn_list,1),length(h_mean_par_period_pre));
tsamplUTC = nan(size(stn_list,1),1);
real_noonLOC = nan(size(stn_list,1),1);
OUT1 = nan(size(stn_list,1),7);

for id = 1:size(stn_list,1)
    
    if icecamp
        tmp = samples(samples(:,1)==year & samples(:,2)==stn_list(id,6),:);
    else
        if size(samples,1) == size(stn_list,1), tmp = samples(id,:);
        else error('Number of casts in samples and stn_list must match in the Amundsen data set')
        end
    end
    UTC_sample = nanmean(tmp(:,5));
    if isnan(UTC_sample), UTC_sample = nanmean(samples(:,5)); end
    
    if ~isempty(tmp) && ~isnan(UTC_sample)
        sl = stn_list(id,:);
        
        % Output for ice camp
        OUT1(id,:) = [sl(5) sl(3) sl(4) sl(6) UTC_sample/24 sl(1) sl(2)];
        
        tsamplUTC(id) = datenum(sl(5),sl(3),sl(4),UTC_sample,0,0);
        real_diff = -sl(2)/15;
        real_noonLOC(id) = 0.5 + (real_diff-UTC_diff)/24;
        
        % In situ PAR data
        PARTMP = mean_par_period(par_ts, tsamplUTC(id), tsearch_h, pXh, nanmean(tmp(:,7)), UTC_diff, tm_insitu);
        PAR0_insitu(id,:) = (cell2mat(struct2cell(PARTMP)))';
        
        % Satellite SBDART-derived data
        EdSpec = [];
        mtimeUTC = [];
        for k = 1:length(backsearch)
            tmpEdSpec = dlmread(sprintf('%s/Ed0_%4i_%03i_%0.3f_%0.3f.txt',datapath,sl(5),sl(6)+backsearch(k),sl(1),sl(2)));
            EdSpec = [EdSpec; tmpEdSpec'];
            mtimeUTC = [mtimeUTC; datenum(sl(5),0,sl(6)+backsearch(k),hours,0,0)];
        end
        par_ts_sbdart.data = mean(EdSpec(:,wlint),2)*wlrange;
        par_ts_sbdart.mtimeUTC = mtimeUTC;
        par_ts_sbdart.mtimeLOC = mtimeUTC - UTC_diff/24;
        PARTMP = mean_par_period(par_ts_sbdart, tsamplUTC(id), tsearch_h, pXh, nanmean(tmp(:,7)), UTC_diff, tm_sbdart);
        PAR0_sbdart(id,:) = (cell2mat(struct2cell(PARTMP)))';
        
    end
end

tsamLOC = tsamplUTC - UTC_diff;
real_noonUTC = real_noonLOC + UTC_diff;

%% Correct bias

select4correction = (...
    {'PARsam' 'PARsam_interp' 'PARp1h' 'PARp3h' 'PAR3hSAT'...
    'PARp24h' 'PARp48h' 'PARp24hmax' 'PARdayLOC' 'PARnoonLOC'...
    'PARnoon1hLOC' 'PARnoon3hSAT' 'PARdayLOCmax' 'PARdayUTC'})';
if correct_bias
    for j = 1:length(select4correction)
        k = find(strcmp(select4correction{j},h_mean_par_period_pre));
        if sum(~isnan(PAR0_insitu(:,k)) & ~isnan(PAR0_sbdart(:,k)))~=0
            % Extract stats
            s = f_skill_stats(PAR0_insitu(:,k),PAR0_sbdart(:,k),'lin','stats4diagramsOFF');
            PAR0_sbdart(:,k) = PAR0_sbdart(:,k)*(1 - s.bias/nanmean(PAR0_sbdart(:,k)));
        end
    end
end


%% Select variables for export and corresponding headers

% h_mean_par_period_output = (...
%     {'PARsam' 'PARsam_interp' 'PARp1h' 'Np1h' 'PARp3h' 'Np3h'...
%     'PARp24h' 'Np24h' 'PARp48h' 'Np48h' 'PARdayLOC'  'NdayLOC' 'PARnoon1hLOC' 'Nnoon1hLOC'})';
iselect = [1:6 9:12 14 15 17 18];
h_mean_par_period_output = h_mean_par_period_pre(iselect);
PAR0_insitu_OUT = PAR0_insitu(:,iselect);
PAR0_sbdart_OUT = PAR0_sbdart(:,iselect);

if icecamp
    rmrows = isnan(OUT1(:,1));
    PAR0_insitu_OUT(rmrows,:) = [];
    PAR0_sbdart_OUT(rmrows,:) = [];
    OUT1(rmrows,:) = [];
end


%% Filter out measurements with insufficient data

ccols = [3 5 7 9 11 13];
for j = 1:length(ccols)
    nmax = nanmax(PAR0_insitu_OUT(:,ccols(j)+1));
    rmdata = PAR0_insitu_OUT(:,ccols(j)+1) < 0.8*nmax;
    PAR0_insitu_OUT(rmdata,ccols(j)) = nan;
end


%% Extract stats for table with selected PAR0 metrics
statcols = [5 7 9];
STATSOUT = nan(7,length(statcols));
for j = 1:length(statcols)
    s = f_skill_stats(PAR0_insitu_OUT(:,statcols(j)),PAR0_sbdart_OUT(:,statcols(j)),'lin','stats4diagramsOFF');
    tmp = [s.r s.rms s.mape 100*s.bias/nanmean(PAR0_sbdart_OUT(:,statcols(j))) s.slope s.slopeMA s.N];
    STATSOUT(:,j) = tmp';
end

%% Plot for final check

if icecamp
    xplot = OUT1(:,4);
    xl = 'DOY';
else
    load ~/Desktop/GreenEdge/Irradiance/CASTS_ICE_matched_25-Apr-2017.mat
    xplot = OUT(:,8);
    xl = 'CAST#';
end
ymin = nanmin(nanmin([PAR0_insitu_OUT(:,[1:3 5 7 9 11 13]) PAR0_sbdart_OUT(:,[1:3 5 7 9 11 13])]));
ymax = nanmax(nanmax([PAR0_insitu_OUT(:,[1:3 5 7 9 11 13]) PAR0_sbdart_OUT(:,[1:3 5 7 9 11 13])]));
figure(919)
subplot(121)
plot(xplot,PAR0_insitu_OUT(:,[1:3 5 7 9 11 13]),'-')
axis([nanmin(xplot)-2 nanmax(xplot)+2 ymin ymax])
xlabel(xl), ylabel('PAR (µE m^{-2} s^{-1}')
legend(h_mean_par_period_output([1:3 5 7 9 11 13]),'interpreter','none')
subplot(122)
plot(xplot,PAR0_sbdart_OUT(:,[1:3 5 7 9 11 13]),'-')
axis([nanmin(xplot)-2 nanmax(xplot)+2 ymin ymax])
xlabel(xl)

%% Format output and export as tab-separated values

if icecamp
    h_out1 = {'yearUTC';'monthUTC';'dayUTC';'doyUTC';'timeUTC';'lat';'lon'};
    OUT_IS = [OUT1 PAR0_insitu_OUT];
    OUT_SBDART = [OUT1 PAR0_sbdart_OUT];
    OUT_IS(isnan(OUT_IS(:,1)),:) = [];
    h = [h_out1; h_mean_par_period_output];
    ICorAM = 'IC';
else
    load ~/Desktop/GreenEdge/Irradiance/CASTS_ICE_matched_25-Apr-2017.mat
    ICE = OUT; ICE(:,10) = ICE(:,10)/10;
    h_ice = {'yearUTC';'monthUTC';'dayUTC';'doyUTC';'timeUTC';'lat';'lon';'cast';'z_bottom';'SIC_visual';'IceMatchDist_km';'SIC_minus2d';'SIC_minus1d';'SIC_day'};
    h = [h_ice; h_mean_par_period_output];
    OUT_IS = [ICE PAR0_insitu_OUT];
    OUT_SBDART = [ICE PAR0_sbdart_OUT];
    ICorAM = 'AM';
end

fid1 = fopen(sprintf('%s/PAR0_insitu_dailyMetrics.GE%s.%s.txt',outpath,num2str(year),ICorAM), 'w');
fid2 = fopen(sprintf('%s/PAR0_sbdart_dailyMetrics.GE%s.%s.txt',outpath,num2str(year),ICorAM), 'w');
for j = 1:length(h)
    if j < length(h)
        fprintf(fid1,'%s\t',h{j});
        fprintf(fid2,'%s\t',h{j});
    else
        fprintf(fid1,'%s\n',h{j});
        fprintf(fid2,'%s\n',h{j});
    end
end

if strcmp(specrange,'PAR')
    dlmwrite(sprintf('%s/PAR0_insitu_dailyMetrics.GE%s.%s.txt',outpath,num2str(year),ICorAM),OUT_IS,'delimiter','\t','-append')
    dlmwrite(sprintf('%s/%s0_sbdart_dailyMetrics.GE%s.%s.txt',outpath,specrange,num2str(year),ICorAM),OUT_SBDART,'delimiter','\t','-append')
else
    outpath = '/Users/martigalitapias/Desktop/GreenEdge/Irradiance/ED0_DAILY_METRICS';
    dlmwrite(sprintf('%s/%s0_sbdart_dailyMetrics.GE%s.%s.txt',outpath,specrange,num2str(year),ICorAM),OUT_SBDART,'delimiter','\t','-append')
end
