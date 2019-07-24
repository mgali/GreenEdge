% Test mean_par_period.m on 2016 ice camp PAR data, plot
clc, clear

load ~/Desktop/GreenEdge/Irradiance/PAR_insitu/PAR_Data_GE2016-ICECAMP_cyber.mat
% sample_list = dlmread('input.noclim.GE2016.txt');
sample_list = dlmread('input.noclim.IceCamp20152016.txt');

% General definitions
tsearch_h = 3.1;
pXh = 48;
UTC_diff = 4;
id = 234;
tsamplUTC = datenum(sample_list(id,5),sample_list(id,3),sample_list(id,4),16,20,0);
lon = -63.78953;

%% Function to debug

% GENERAL
its = 1:size(par_ts.data,1); % index list of PAR time series
tsearch = tsearch_h/24; % search period in days (Matlab time is in days)

% SAMPLING TIME AND BACKWARDS-INTEGRATED DATA
tdif = par_ts.mtimeUTC - tsamplUTC; % diff between time series and sampling
iref_all = its(abs(tdif) <= tsearch); % indexs of measurements within tsearch

if iref_all % test that matches were found
    
    % 1. PAR at sampling time interpolated between closest measurements
    % Find closest measurement within ±3 hours (this interval was selected
    % because output of SBDART RT code has 3h resolution). Get instantaneous
    % PAR estimate from linear interpolation of two closest measurements
    
    irefminus = nanmax(iref_all(tdif(iref_all) <= 0));
    irefplus = nanmin(iref_all(tdif(iref_all) > 0));
    
    % If multiple matches found, select the closest, starting with time
    % stamps prior to sampling time since time corresponds to start of PAR integration
    if irefminus
        iref = irefminus;
    elseif irefplus
        iref = irefplus;
    end
    par0plus.tsam = par_ts.data(iref);
    
    % Interpolate between two measurements if possible
    if irefminus && irefplus
        ii = [irefminus irefplus];
        if sum(~isnan(ii)) == 2
            par0plus.tsami = interp1(par_ts.mtimeUTC(ii),par_ts.data(ii),tsamplUTC);
        end
    end
    
    % 2.1 Mean PAR during 24 hours prior to sampling, corresponding data counts
    p24h = 24;
    i24h = (par_ts.mtimeUTC <= par_ts.mtimeUTC(iref)) & (par_ts.mtimeUTC > par_ts.mtimeUTC(iref) - p24h/24);
    par0plus.p24h = nanmean(par_ts.data(i24h));
    par0plus.Np24h = sum(~isnan(par_ts.data(i24h)));
    
    % 2.2 Mean PAR during X hours prior to sampling, corresponding data counts
    if ~isnan(pXh) && ~isempty(pXh)
        iXh = (par_ts.mtimeUTC <= par_ts.mtimeUTC(iref)) & (par_ts.mtimeUTC > par_ts.mtimeUTC(iref) - pXh/24);
        par0plus.pXh = nanmean(par_ts.data(iXh));
        par0plus.NpXh = sum(~isnan(par_ts.data(iXh)));
    else
        par0plus.pXh = nan;
        par0plus.NpXh = nan;
    end
    
else
    
    par0plus.tsam = nan;
    par0plus.tsami = nan;
    par0plus.p24h = nan;
    par0plus.Np24h = nan;
    par0plus.pXh = nan;
    par0plus.NpXh = nan;
    
end

% NATURAL DAY DATA
t_vec = datevec(par_ts.mtimeLOC);
tsamplLOCAL = tsamplUTC - UTC_diff/24;
ts_vec = datevec(tsamplLOCAL);
iday = (t_vec(:,1) == ts_vec(1)) & (t_vec(:,2) == ts_vec(2)) & (t_vec(:,3) == ts_vec(3));

if sum(iday)
    
    % 3. Mean daily PAR during the LOCAL day, data counts
    data_iday = par_ts.data(iday);
    par0plus.day = nanmean(data_iday);
    par0plus.Nday = sum(~isnan(data_iday));
    
    % 4. Mean PAR at solar noon on sampling day
    % Implement interpolations based on interp_PAR_tests.m
    
    % 4.1 Take value from closest measurement period (backwards temporal search
    % since we look at integration that started before noon time, not after)
    real_diff = -lon/15;
    UTC_noon = 12 + real_diff;
    mtimeUTC_iday = par_ts.mtimeUTC(iday);
    mtimeUTC_iday_vec = datevec(mtimeUTC_iday);
    hourdiff = mtimeUTC_iday_vec(:,4) + mtimeUTC_iday_vec(:,5)/60 - UTC_noon;
    mhd = nanmax(hourdiff(hourdiff<=0));
    if ~isnan(mhd) && (abs(mhd)/24 < tsearch)
        par0plus.noon = data_iday(hourdiff==mhd);
    else
        par0plus.noon = nan;
    end
    
else
    
    par0plus.day = nan;
    par0plus.Nday = nan;
    par0plus.noon = nan;
    
end
