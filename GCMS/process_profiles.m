% PROCESS DMS PROFILES FROM GREEN EDGE
clear
close all
% data = xlsread('~/Desktop/GreenEdge/GCMS/gcms_greenedge_Marti.xlsx','Niskin_underway','A2:AD139'); % 30 columns
data_dmspt = xlsread('~/Desktop/GreenEdge/GCMS/gcms_greenedge_Marti.xlsx','Niskin_underway','A2:AL140');
data = data_dmspt(:,[1:30 end]);

% Sort
[data, sortindexs] = sortrows(data,[1 2 3]);

% Date used to extract daily data in loop
date = datenum(data(:,1:3));
diffdate = diff(date);
dateax = [date(diffdate~=0); date(end)];
dateaxstr = cellstr(datestr(dateax));

% Prepare time-dependent calibration: daily interpolation of [doy slope intercept]
cal_td_pre = [
    2016 0 176 0 0 0       0.7193	12.1171
    2016 0 180 0 0 0       0.968     13.262
    2016 0 183.01 0 0 0	  1.0703	13.568
    2016 0 183.02 0 0 0	 0.9887	13.266
    2016 0 187 0 0 0   	0.9955	13.464
    2016 0 192 0 0 0       0.9577	13.006];
cal_td = interp1(datenum(cal_td_pre(:,1:6)),cal_td_pre(:,7:8),dateax);

% Define opt1 calibration (expected concentration in 20 mL sample w/r/to peak area)
% cal_opt1 = [0.968 13.26]; % doy 180 calibration
cal_opt1 = [0.9775 13.2495]; % mean of 4 more similar calibrations, nanmean(cal_td_pre([2 4:6],end-1:end))

% Define opt2 calibration
cal_opt2 = [0.8402 2.31]; % pmol S w/r/to peak area, no data excluded
% cal_opt2 = [0.8527 2.509]; % pmol S w/r/to peak area, 24/06/2016 calibration excluded

% Define baseline correction for residual DMS 62 that remains in the system despite MQ rinsing (calculated in spreadsheet).
% This amount of S Will be removed AFTER applying calibration to 62 peak area.
% Otherwise we would be applying calibration to too small peak areas with larger associated uncertainty/
% /////////// REFINE, MAKE TIME-DEPENDENT ///////////
corr_pa_fixed = 10000; % minimal 62 peak area found in injections where 62 was not expected (24000?)
% ///////////////////////////////////////////////////

% Define correction factor for samples stored frozen (transect 400)
corr_fs = 1; % SET TO 1: Not needed when applying internal standard correction on a per sample basis!


%% LOOP FOR EXTRACTING DAILY PROFILE DATA, APPLYING CORRECTIONS/CALIBRATIONS
% AND EXPORTING CORRECTED/CALIBRATED DATA

DMSOUT = [];
headout = ({'dms_td','dms_opt1','dms_opt2','dms_cal2point',...
    'dms_consens','e_dms_consens','cv_dms_consens','cf68_td','cf68_opt1','cf68_opt2',...
    'dms_td_cf68','dms_opt1_cf68','dms_opt2_cf68','dms_cal2point_cf68',...
    'dms_consens_cf68','e_dms_consens_cf68','cv_dms_consens_cf68','dmspt'})';

for j = 1:length(dateax)
    
    TMP = data(date == dateax(j),:);
    OUT = nan(size(TMP,1),9);
    
    % Assign variable names: do it all inside loop!
    %year = TMP(:,1);
    %month = TMP(:,2);
    %day = TMP(:,3);
    %time_dec = TMP(:,4);
    %station = TMP(:,6);
    %z = TMP(:,13);
    injvol = TMP(:,14); % injection volume in mL
    %prescreen = TMP(:,15); % prescreening applied (100 vs 5 µm mesh)
    pa62 = TMP(:,16); % peak areas
    pa65 = TMP(:,17);
    pa68 = TMP(:,18);
    ec65 = TMP(:,25); % expected concentrations in nmol/L
    ec68 = TMP(:,26);
    correct_frozen = nan(size(TMP,1),1);
    correct_frozen(TMP(:,29)==0) = 1;
    correct_frozen(TMP(:,29)~=0) = corr_fs;
    
    % ///// Baseline DMS correction ///// 
%     corr_pa = corr_pa_fixed*ones(size(TMP,1),1); % 1. Fixed corr_pa
%     corr_pa = 0.1*nanmedian(pa68)*ones(size(TMP,1),1); % 2. Fixed daily corr_pa
    corr_pa = 0.1*pa68; % 3. Per sample corr_pa
    
    % Calculate DMS with time-dependent calibration
    cal_day = cal_td(j,:);
    correct_td = 1e9*(20./injvol).*10.^(cal_day(1)*log10(corr_pa)-cal_day(2)).*correct_frozen;
    OUT(:,1) = (1e9*(20./injvol).*10.^(cal_day(1)*log10(pa62)-cal_day(2)) - correct_td).*correct_frozen;
    meas68_td = 1e9*(20./injvol).*10.^(cal_day(1)*log10(pa68)-cal_day(2)).*correct_frozen;
    OUT(:,8) = ec68./meas68_td;
    
    % Calculate DMS with unique optimized calibration subjectively chosen
    correct_opt1 = 1e9*(20./injvol).*10.^(cal_opt1(1)*log10(corr_pa)-cal_opt1(2)).*correct_frozen;
    OUT(:,2) = (1e9*(20./injvol).*10.^(cal_opt1(1)*log10(pa62)-cal_opt1(2)) - correct_opt1).*correct_frozen;
    meas68_opt1 = 1e9*(20./injvol).*10.^(cal_opt1(1)*log10(pa68)-cal_opt1(2)).*correct_frozen;
    OUT(:,9) = ec68./meas68_opt1;
    
    % Calculate DMS with unique optimized calibration objectively chosen
    correct_opt2 = (10.^(cal_opt2(1)*log10(corr_pa)-cal_opt2(2))./injvol).*correct_frozen;
    OUT(:,3) = (10.^(cal_opt2(1)*log10(pa62)-cal_opt2(2))./injvol - correct_opt2).*correct_frozen;
    meas68_opt2 = (10.^(cal_opt2(1)*log10(pa68)-cal_opt2(2))./injvol).*correct_frozen;
    OUT(:,10) = ec68./meas68_opt2;
    
    % Calculate DMS with sample-wise interpolation between 2 internal
    % standards... less robust than other 3, increases the CV.
    OUT(:,4) = cal_2point(TMP,correct_frozen,nanmedian(corr_pa));
    
    % Set negative values to zero
    OUT(OUT<0) = 0;
    
    % Append
    DMSOUT = [DMSOUT; OUT];
    
end


%% Apply correction factor based on DMS 68 and calculate std and CV

col2corr = 11:13;
col2notcorr = 14;
DMSOUT(:,col2corr) = DMSOUT(:,col2corr-10).*DMSOUT(:,col2corr-3);
DMSOUT(:,col2notcorr) = DMSOUT(:,col2notcorr-10);

% Stats for uncorrected DMS
DMSOUT(:,5) = (nanmean((DMSOUT(:,1:2))'))';
%DMSOUT(:,6) = (nanstd((DMSOUT(:,1:3))'))';
DMSOUT(:,6) = (nanstd((DMSOUT(:,1:2))'))'/sqrt(2);
DMSOUT(:,7) = DMSOUT(:,6)./DMSOUT(:,5);

% Stats for 68-corrected DMS
DMSOUT(:,15) = (nanmean((DMSOUT(:,11:12))'))';
%DMSOUT(:,16) = (nanstd((DMSOUT(:,11:13))'))';
DMSOUT(:,16) = (nanstd((DMSOUT(:,11:12))'))'/sqrt(2);
DMSOUT(:,17) = DMSOUT(:,16)./DMSOUT(:,15);

% DMSPt
DMSOUT(:,18) = data(:,end);

%% Compare 3 good calibrations after 68 correction and choose the 2 that are more similar

% figure(), set(gcf,'units','centimeters','position',[5 10 55 15])
% subplot(131)
% plot([0.1 80],[0.1 80],'-r','linewidth',2), axis([0.1 80 0.1 80]), axis square, hold on
% plot(DMSOUT(:,11),DMSOUT(:,12),'ok')
% set(gca,'xscale','log','yscale','log'), legend('1:1',sprintf('r = %0.3f',corr(DMSOUT(:,11),DMSOUT(:,12),'rows','pairwise')),'location','northwest')
% subplot(132)
% plot([0.1 80],[0.1 80],'-r','linewidth',2), axis([0.1 80 0.1 80]), axis square, hold on
% plot(DMSOUT(:,11),DMSOUT(:,13),'ok')
% set(gca,'xscale','log','yscale','log'), legend('1:1',sprintf('r = %0.3f',corr(DMSOUT(:,11),DMSOUT(:,13),'rows','pairwise')),'location','northwest')
% subplot(133)
% plot([0.1 80],[0.1 80],'-r','linewidth',2), axis([0.1 80 0.1 80]), axis square, hold on
% plot(DMSOUT(:,12),DMSOUT(:,13),'ok')
% set(gca,'xscale','log','yscale','log'), legend('1:1',sprintf('r = %0.3f',corr(DMSOUT(:,12),DMSOUT(:,13),'rows','pairwise')),'location','northwest')

%% Compare CV for the 2-3 calibrations before and after DMS 68-based cf

compare_cv = [nanmean(DMSOUT(:,7)) nanmean(DMSOUT(:,17));
    nanmin(DMSOUT(:,7)) nanmin(DMSOUT(:,17));
    quantile(DMSOUT(:,7),.05) quantile(DMSOUT(:,17),.05);
    quantile(DMSOUT(:,7),.25) quantile(DMSOUT(:,17),.25);
    nanmedian(DMSOUT(:,7)) nanmedian(DMSOUT(:,17));
    quantile(DMSOUT(:,7),.75) quantile(DMSOUT(:,17),.75);
    quantile(DMSOUT(:,7),.95) quantile(DMSOUT(:,17),.95);
    nanmax(DMSOUT(:,7)) nanmax(DMSOUT(:,17))];

compare_std = [nanmean(DMSOUT(:,6)) nanmean(DMSOUT(:,16));
    nanmin(DMSOUT(:,6)) nanmin(DMSOUT(:,16));
    quantile(DMSOUT(:,6),.05) quantile(DMSOUT(:,16),.05);
    quantile(DMSOUT(:,6),.25) quantile(DMSOUT(:,16),.25);
    nanmedian(DMSOUT(:,6)) nanmedian(DMSOUT(:,16));
    quantile(DMSOUT(:,6),.75) quantile(DMSOUT(:,16),.75);
    quantile(DMSOUT(:,6),.95) quantile(DMSOUT(:,16),.95);
    nanmax(DMSOUT(:,6)) nanmax(DMSOUT(:,16))];

% Plot correction factors
figure(),plot(date,DMSOUT(:,8:10),'o'), title('68 correction factor vs date')
figure(),plot(DMSOUT(:,5),DMSOUT(:,8:10),'o'), title('68 correction factor vs uncorr mean DMS')
figure(),plot(DMSOUT(:,15),DMSOUT(:,8:10),'o'), title('68 correction factor vs corr mean DMS')

% Compare corrected data
% Log version
figure(), axis square, hold on
plot(DMSOUT(:,15),zeros(size(DMSOUT(:,15))),'-r','linewidth',2)
errorbar(DMSOUT(:,15),log10(DMSOUT(:,5)./DMSOUT(:,15)),log10(DMSOUT(:,6)./DMSOUT(:,15)),'.','markersize',20)
% Linear version
figure(), axis square, hold on
plot(DMSOUT(:,15),ones(size(DMSOUT(:,15))),'-r','linewidth',2)
errorbar(DMSOUT(:,15),DMSOUT(:,5)./DMSOUT(:,15),DMSOUT(:,6)./DMSOUT(:,15),'.','markersize',20)

%% Save data

% Data to put back in spreadsheet
PUTBACK = nan(size(DMSOUT,1),11);
iputback = [1 2 11 12 15 16];
PUTBACK(sortindexs,:) = [sortindexs data(:,[1:3 6]) DMSOUT(:,iputback)];

h_putback = headout(iputback);

% Data in *.mat format
save('gcms_greenedge_proc.mat','DMSOUT','headout','data','date','dateax','sortindexs')


