% Test mean_par_period.m on 2016 ice camp PAR data, plot
clc, clear


%% ---------------------------------- EDIT --------------------------------

% In situ dataset: 0 is ice camp, 1 is Amundsen
% icecamp = 0; year = 2016;
% icecamp = 1; year = 2015;
icecamp = 1; year = 2016;

% mean_par_period.m settings
tsearch_h = 3.1;
pXh = 48;
UTC_diff = 4;

% SBDART settings
hperiod = 1;
specres = 5;
wl = 290:specres:700;
wlpar = wl>=400 & wl<=700;
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
        corrfact = 1;
    elseif year == 2016
        load ~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2016-ICECAMP_cyber.mat
        samples = samplesIC2016;
        header_samples = headerIC2016;
        tm_insitu = 1;
        corrfact = 1;
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
    corrfact = 1.00; % 1.06
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


%% Preallocate and fill in situ and SBDART PAR data

% Assign variables
h_mean_par_period_output = (...
    {'PARsam' 'PARsam_interp' 'PARp1h' 'Np1h' 'PARp3h' 'Np3h' 'PAR3hSAT' 'N3hSAT'...
    'PARp24h' 'Np24h' 'PARp48h' 'Np48h' 'PARp24hmax' 'PARdayLOC' 'NdayLOC' 'PARnoonLOC'...
    'PARnoon1hLOC' 'Nnoon1hLOC' 'PARnoon3hSAT' 'Nnoon3hSAT' 'PARdayLOCmax' 'PARdayUTC' 'NdayUTC'})';

PAR0_insitu = nan(size(stn_list,1),length(h_mean_par_period_output));
PAR0_sbdart = nan(size(stn_list,1),length(h_mean_par_period_output));
tsamplUTC = nan(size(stn_list,1),1);
real_noonLOC = nan(size(stn_list,1),1);

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
        par_ts_sbdart.data = mean(EdSpec(:,wlpar),2)*300*corrfact;
        par_ts_sbdart.mtimeUTC = mtimeUTC;
        par_ts_sbdart.mtimeLOC = mtimeUTC - UTC_diff/24;
        PARTMP = mean_par_period(par_ts_sbdart, tsamplUTC(id), tsearch_h, pXh, nanmean(tmp(:,7)), UTC_diff, tm_sbdart);
        PAR0_sbdart(id,:) = (cell2mat(struct2cell(PARTMP)))';
        
    end
end

tsamLOC = tsamplUTC - UTC_diff;
real_noonUTC = real_noonLOC + UTC_diff;

%% Correct bias

% select4correction = (...
%     {'PARsam' 'PARsam_interp' 'PARp1h' 'PARp3h' 'PAR3hSAT'...
%     'PARp24h' 'PARp48h' 'PARp24hmax' 'PARdayLOC' 'PARnoonLOC'...
%     'PARnoon1hLOC' 'PARnoon3hSAT' 'PARdayLOCmax' 'PARdayUTC'})';
select4correction = (...
    {'PARp1h' 'PARp3h' 'PARp24h' 'PARp48h' 'PARnoon1hLOC' 'PARdayLOC'})';
select4correction = (...
    {'PARp3h'})';

if correct_bias
    for j = 1:length(select4correction)
        k = find(strcmp(select4correction{j},h_mean_par_period_output));
        if sum(~isnan(PAR0_insitu(:,k)) & ~isnan(PAR0_sbdart(:,k)))~=0
            % Extract stats
            s = f_skill_stats(PAR0_insitu(:,k),PAR0_sbdart(:,k),'lin','stats4diagramsOFF');
            PAR0_sbdart(:,k) = PAR0_sbdart(:,k)*(1 - s.bias/nanmean(PAR0_sbdart(:,k)));
        end
    end
end

%% Filter out measurements with insufficient data

ccols = [ 3 5 7 9 11 14 17 19 22];
for j = 1:length(ccols)
    nmax = nanmax(PAR0_insitu(:,ccols(j)+1));
    rmdata = PAR0_insitu(:,ccols(j)+1) < 0.8*nmax;
    PAR0_insitu(rmdata,ccols(j)) = nan;
end


%% --------- Plot original time series along with extracted values --------

% Clean matrices
toremove = sum(isnan(PAR0_sbdart),2)==size(PAR0_sbdart,2) | sum(isnan(PAR0_insitu),2)==size(PAR0_insitu,2);
PAR0_insitu(toremove,:) = [];
PAR0_sbdart(toremove,:) = [];
stn_list(toremove,:) = [];

% Select variables to plot
% vnames = ({'PARsam' 'PARp1h' 'PARp3h'})'; vcolors = [0 0 1; 0 .6 0; 1 0 0];
% vnames = ({'PARp24h' 'PARp48h' 'PARdayLOC'})'; vcolors = [0 0 1; 0 .6 0; 1 0 0];
vnames = ({'PARp24h' 'PARp48h'})'; vcolors = [0 0 1; 0 .6 0; 1 0 0];
% vnames = ({'PARnoonLOC' 'PARnoon1hLOC' 'PARdayLOCmax'})'; vcolors = [0 0 1; 0 .6 0; 1 0 0];
% vnames = ({'PARp3h' 'PARp24h'})'; vcolors = [0 0 1; 0 .6 0; 1 0 0];
sc = 'lin'; % 'lin' or 'log'
X = nan(size(stn_list,1),length(vnames));
Y = nan(size(X));

if icecamp
    % Scatterplot, first impression, compare different periods
    figure(111), clf, hold on
    title('PAR (µmol quanta m^{-2} s^{-1})','fontsize',20)
    set(gcf,'units','centimeters','position',[3 2 20 20])
    
    for j = 1:length(vnames)
        X(:,j) = PAR0_insitu(:,strcmp(vnames{j},h_mean_par_period_output));
        Y(:,j) = PAR0_sbdart(:,strcmp(vnames{j},h_mean_par_period_output));
        
        figure(111)
        scatter(X(:,j),Y(:,j),50,vcolors(j,:),'linewidth',2)
    end
    
    figure(111),
    legend(vnames,'location','southeast')
    Xc = X(:); Yc = Y(:);
    if strcmp(sc,'lin'), xymin = nanmin([Xc; Yc]);
    elseif strcmp(sc,'log'), xymin = nanmin([Xc(Xc>0 & Yc>0); Yc(Xc>0 & Yc>0)]);
    end
    xymax = nanmax([X(:); Y(:)]);
    axis([xymin xymax xymin xymax])
    plot([xymin xymax],[xymin xymax],'-k')
    plot([xymin xymax],1.5*([xymin xymax]),'--k')
    plot([xymin xymax],(1/1.5)*([xymin xymax]),'--k')
    plot([xymin xymax],2*([xymin xymax]),'-.k')
    plot([xymin xymax],.5*([xymin xymax]),'-.k')
    plot([xymin xymax],3*([xymin xymax]),':k')
    plot([xymin xymax],(1/3)*([xymin xymax]),':k')
    xlabel('In situ','fontsize',20)
    ylabel('SBDART-Aqua','fontsize',20)
    set(gca,'fontsize',16,'xscale',sc,'yscale',sc)
    % Print stats on graph
    xt = xymax-xymin;
    yt = xymax-xymin;
    for j = 1:length(vnames)
        s = f_skill_stats(X(:,j),Y(:,j),'off','off');
        if ~isempty(s)
            text(xymin+xt*(0.02+0.23*(j-1)),xymax-yt*0.03,sprintf('r = %0.2f',s.r),'fontsize',14,'color',vcolors(j,:))
            text(xymin+xt*(0.02+0.23*(j-1)),xymax-yt*0.08,sprintf('RMSE = %0.0f',s.rms),'fontsize',14,'color',vcolors(j,:))
            text(xymin+xt*(0.02+0.23*(j-1)),xymax-yt*0.13,sprintf('Bias = %0.0f%%',100*s.bias/nanmean(X(:,j))),'fontsize',14,'color',vcolors(j,:))
            text(xymin+xt*(0.02+0.23*(j-1)),xymax-yt*0.18,sprintf('MAPE = %0.1f%%',s.mape),'fontsize',14,'color',vcolors(j,:))
            text(xymin+xt*(0.02+0.23*(j-1)),xymax-yt*0.23,sprintf('Slope = %0.2f',s.slope),'fontsize',14,'color',vcolors(j,:))
            text(xymin+xt*(0.02+0.23*(j-1)),xymax-yt*0.28,sprintf('SlopeMA = %0.2f',s.slopeMA),'fontsize',14,'color',vcolors(j,:))
            text(xymin+xt*(0.02+0.23*(j-1)),xymax-yt*0.33,sprintf('N = %i',s.N),'fontsize',14,'color',vcolors(j,:))
        end
    end
    
    % Ratio, for icecamp
    figure(112), clf, hold on
    set(gcf,'units','centimeters','position',[3 25 30 20])
    xax = stn_list(:,6);
    bar(xax,Y./X)
    xlim([nanmin(xax)-2 nanmax(xax)+2])
    ylim([0 2])
    legend(vnames,'location','northwest','fontsize',14)
    plot([nanmin(xax)-2 nanmax(xax)+2],[1 1],'-k')
    xlabel(sprintf('Day of year %4i',year),'fontsize',20)
    ylabel('Ratio SBDART-Aqua / In situ','fontsize',20)
    set(gca,'fontsize',16)
end

%% GE2016 cruise only, compare deviations between SBDART and data as a function of ice cover

if ~icecamp
    
    load ~/Desktop/GreenEdge/Irradiance/CASTS_ICE_matched_24-Mar-2017.mat
    iceday = OUT(:,strcmp('SIC_day',header_out));
    icecut = 0.75;
    
    vnames = select4correction;
    %     vnames = ({'PARsam'})';
    %     vnames = ({'PARp1h'})';
    %     vnames = ({'PARp3h'})';
    %     vnames = ({'PARp24h'})'; vcolors = [0 0 1; 0 .6 0; 1 0 0]; iceday = OUT(:,strcmp('SIC_day',header_out));
    %     vnames = ({'PARp48h'})'; vcolors = [0 0 1; 0 .6 0; 1 0 0]; iceday = nanmean(OUT(:,(end-1):end),2);
    sc = 'lin'; % 'lin' or 'log'
    
    for j = 1:length(vnames)
        
        X = PAR0_insitu(:,strcmp(vnames{j},h_mean_par_period_output));
        Y = PAR0_sbdart(:,strcmp(vnames{j},h_mean_par_period_output));
        Xc = X(:); Yc = Y(:);
        if strcmp(sc,'lin'), xymin = nanmin([Xc; Yc]);
        elseif strcmp(sc,'log'), xymin = nanmin([Xc(Xc>0 & Yc>0); Yc(Xc>0 & Yc>0)]);
        end
        xymax = nanmax([X(:); Y(:)]);
        
        figure(j*10), clf, hold on
        title('PAR (µmol quanta m^{-2} s^{-1})','fontsize',20)
        set(gcf,'units','centimeters','position',[2 3 50 15])
        
        subplot(131)
        scatter(X(iceday<=icecut),Y(iceday<=icecut),50,'b','linewidth',2), hold on
        scatter(X(iceday>icecut),Y(iceday>icecut),50,[1 .5 0],'linewidth',2)
        
        subplot(132)
        scatter(iceday(iceday<=icecut),Y(iceday<=icecut)./X(iceday<=icecut),50,'b','linewidth',2), hold on
        scatter(iceday(iceday>icecut),Y(iceday>icecut)./X(iceday>icecut),50,[1 .5 0],'linewidth',2)
        
        subplot(133)
        binice = (floor((iceday/2)*10)+1)/5 - 0.1;
        yb2 = Y./X;
        plot([0 6],[1 1],'-k','linewidth',2), hold on
        plot([xymin xymax],[xymin xymax],'-k')
        plot([xymin xymax],1.5*([xymin xymax]),':k')
        plot([xymin xymax],0.5*([xymin xymax]),':k')
        axis([xymin xymax xymin xymax])
        %         axis square
        boxplot(yb2,binice)
        grid on
        set(gca,'xtick',1.5:5.5,'xticklabel',0.2:0.2:1,'xscale',sc,'yscale',sc,...
            'tickdir','out','ticklength',[.02 .02],'fontsize',16)
        xlabel(sprintf('\nSea ice concentration'),'fontsize',18)
        ylabel('SBDART / In situ','fontsize',18)
        ylim([0 2])
        
        subplot(131)
        xlabel('PAR In situ, µE m^{-2} s^{-1}','fontsize',18)
        ylabel('PAR SBDART-Aqua, µE m^{-2} s^{-1}','fontsize',18)
        set(gca,'fontsize',16,'xscale',sc,'yscale',sc,'tickdir','out','ticklength',[.02 .02])
        Xc = X(:); Yc = Y(:);
        if strcmp(sc,'lin'), xymin = nanmin([Xc; Yc]);
        elseif strcmp(sc,'log'), xymin = nanmin([Xc(Xc>0 & Yc>0); Yc(Xc>0 & Yc>0)]);
        end
        xymax = nanmax([X(:); Y(:)]);
        plot([xymin xymax],[xymin xymax],'-k')
        plot([xymin xymax],1.5*([xymin xymax]),':k')
        plot([xymin xymax],0.5*([xymin xymax]),':k')
        axis([xymin xymax xymin xymax]), axis square
        
        % Print stats on graph
        xt = xymax-xymin;
        yt = xymax-xymin;
        s = f_skill_stats(X(iceday<=icecut),Y(iceday<=icecut),'lin','stats4diagramsON');
        if ~isempty(s)
            text(xymin+xt*0.02,xymax-yt*0.05,sprintf('r = %0.2f',s.r),'fontsize',12,'color','b')
            text(xymin+xt*0.02,xymax-yt*0.1,sprintf('RMSE = %0.0f',s.rms),'fontsize',12,'color','b')
            % text(xymin+xt*0.02,xymax-yt*0.15,sprintf('rRMSE = %0.0f%%',100*s.rms_star),'fontsize',12,'color','b')
            text(xymin+xt*0.02,xymax-yt*0.15,sprintf('MAPE = %0.0f',s.mape),'fontsize',12,'color','b')
            text(xymin+xt*0.02,xymax-yt*0.2,sprintf('Bias = %0.0f%%',100*s.bias/nanmean(X(iceday<=icecut))),'fontsize',12,'color','b')
            text(xymin+xt*0.02,xymax-yt*0.25,sprintf('Slope = %0.2f',s.slope),'fontsize',12,'color','b')
            % text(xymin+xt*0.02,xymax-yt*0.35,sprintf('SlopeMA = %0.2f',s.slopeMA),'fontsize',12,'color','b')
            text(xymin+xt*0.02,xymax-yt*0.3,sprintf('N = %i',s.N),'fontsize',12,'color','b')
        end
        s = f_skill_stats(X(iceday>icecut),Y(iceday>icecut),'lin','stats4diagramsON');
        if ~isempty(s)
            text(xymin+xt*0.26,xymax-yt*0.05,sprintf('r = %0.2f',s.r),'fontsize',12,'color',[1 .5 0])
            text(xymin+xt*0.26,xymax-yt*0.1,sprintf('RMSE = %0.0f',s.rms),'fontsize',12,'color',[1 .5 0])
            % text(xymin+xt*0.26,xymax-yt*0.15,sprintf('rRMSE = %0.0f%%',100*s.rms_star),'fontsize',12,'color',[1 .5 0])
            text(xymin+xt*0.26,xymax-yt*0.15,sprintf('MAPE = %0.0f',s.mape),'fontsize',12,'color',[1 .5 0])
            text(xymin+xt*0.26,xymax-yt*0.2,sprintf('Bias = %0.0f%%',100*s.bias/nanmean(X(iceday>icecut))),'fontsize',12,'color',[1 .5 0])
            text(xymin+xt*0.26,xymax-yt*0.25,sprintf('Slope = %0.2f',s.slope),'fontsize',12,'color',[1 .5 0])
            % text(xymin+xt*0.26,xymax-yt*0.35,sprintf('SlopeMA = %0.2f',s.slopeMA),'fontsize',12,'color',[1 .5 0])
            text(xymin+xt*0.26,xymax-yt*0.3,sprintf('N = %i',s.N),'fontsize',12,'color',[1 .5 0])
        end
        
        subplot(132)
        plot([0 1],[1 1],'-k')
        plot([0 1],[1.5 1.5],':k')
        plot([0 1],[.5 .5],':k')
        axis([0 1 0 2]), axis square
        xlabel('Sea ice concentration','fontsize',18)
        ylabel('SBDART / In situ','fontsize',18)
        set(gca,'fontsize',16,'xscale',sc,'yscale',sc,'tickdir','out','ticklength',[.02 .02])
        vv = ({vnames{j}})';
        ll = [vv; 'Idem Ice>0.75'; '1:1'; '±50%'];
        legend(ll,'location','north')
        
    end
    
end

