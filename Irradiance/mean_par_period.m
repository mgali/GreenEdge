% FUNCTION mean_par_period.m
%
% Mart� Gal� T�pias 20170201
%
% INPUTS:
%   par_ts: a PAR time series with matlab-format local and UTC time and PAR data
%   tsamplUTC: the sampling date and time, UTC
%   tsearch: the search time interval
%   nhours: the hours prior to sampling for which the mean PAR is calculated
%   (the prior 24h average is calculated by default)
%   pXh: the hours prior to sampling for which the mean PAR is calculated
%   lon: longitude
%   UTC_diff: difference between UTC and local time, integer hours (e.g. -4)
%   tm:  measurement period in hours. Usually 3 (SBDART), 1 (IceCamp2016),
%   1/60 (Amundsen2016)
%
% Ideally, there should be no gaps in the PAR time series. Nan will be
% omitted
%
% OUTPUTS:
%   par0plus.tsam: PAR closest to sampling time
%   par0plus.tsami: PAR at sampling time interpolated between two closest measurements
%   par0plus.tsam1h: PAR during UTC hour prior to sampling (and counts, N)
%   par0plus.tsam1h: PAR during 3 hours prior to sampling (and counts, N)
%   par0plus.tsam3hSAT: PAR during UTC 3-hour period comprising sampling time (and counts, N)
%   (0-3-6-9-12-15-18-21 as in SBDART)
%   par0plus.p24h: Mean PAR during n hours prior to sampling
%   par0plus.Np24h: Count data used to compute par0plus.pXh
%   par0plus.pXh: Mean PAR during n hours prior to sampling
%   par0plus.NpXh: Count data used to compute par0plus.pXh
%   par0plus.dayLOC: Mean daily PAR irrespective of sampling time, LOCAL day
%   par0plus.Nday: Count data used to compute par0plus.day
%   par0plus.noonLOC: Mean PAR at solar noon on sampling day (different
%       methods might be used, see code interp_PAR_tests.m
%   par0plus.noon1h: Mean 1-h PAR at solar noon (and counts, N)
%   par0plus.noon3h: Mean 3-h PAR at solar noon (0-3-6-9-12-15-18-21 as in
%   SBDART)  (and counts, N)
%
% It is generally assumed that the input time series have been
% quality-controlled, that there are minimal or no gaps, and that they are
% regularly spaced at a known measurement period

function par0plus = mean_par_period(par_ts, tsamplUTC, tsearch_h, pXh, lon, UTC_diff, tm)

if nargin ~= 7
    error('Must supply correct number of arguments')
end

% INITIALIZE ALL OUTPUT WITH NAN
par0plus.tsam = nan;
par0plus.tsami = nan;
par0plus.tsam1h = nan;
par0plus.Ntsam1h = nan;
par0plus.tsam3h = nan;
par0plus.Ntsam3h = nan;
par0plus.tsam3hSAT = nan;
par0plus.Ntsam3hSAT = nan;
par0plus.p24h = nan;
par0plus.Np24h = nan;
par0plus.pXh = nan;
par0plus.NpXh = nan;
par0plus.p24hmax = nan;
par0plus.dayLOC = nan;
par0plus.NdayLOC = nan;
par0plus.noonLOC = nan;
par0plus.noon1hLOC = nan;
par0plus.Nnoon1hLOC = nan;
par0plus.noon3hSATLOC = nan;
par0plus.Nnoon3hSATLOC = nan;
par0plus.dayLOCmax = nan;
par0plus.dayUTC = nan;
par0plus.NdayUTC = nan;
% par0plus.dayLOCmaxUTC = nan; % if eventually no data found, nanmax([])=[],
% whereas nanmean([]) = nan; this can cause problems

% GENERAL
its = 1:size(par_ts.data,1); % index list of PAR time series
tsearch = tsearch_h/24; % search period in days (Matlab time is in days)

% SAMPLING TIME AND BACKWARDS-INTEGRATED DATA
tdif = par_ts.mtimeUTC - tsamplUTC; % diff between time series and sampling
iref_all = its(abs(tdif) <= tsearch); % indexs of measurements within tsearch
iref_all(isnan(iref_all)) = [];

if iref_all % test that matches were found
    
    % 1. PAR at sampling time interpolated between closest measurements
    % Find closest measurement within �3 hours (this interval was selected
    % because output of SBDART RT code has 3h resolution). Get instantaneous
    % PAR estimate from linear interpolation of two closest measurements
    
    irefminus = nanmax(iref_all(tdif(iref_all) <= 0));
    irefplus = nanmin(iref_all(tdif(iref_all) > 0));
    
    % If multiple matches found, select the closest, starting with time
    % stamps prior to sampling time since time corresponds to start of PAR integration
    if ~isempty(irefminus)
        iref = irefminus;
    elseif ~isempty(irefplus)
        iref = irefplus;
    end
    par0plus.tsam = par_ts.data(iref);
    par0plus.tsam3hSAT = par0plus.tsam; % assign tsam value temporarily
    par0plus.tsam1h = par0plus.tsam; % assign tsam value temporarily
    
    % Interpolate between two measurements if possible
    if ~isempty(irefminus) && ~isempty(irefplus)
        ii = [irefminus irefplus];
        par0plus.tsami = interp1(par_ts.mtimeUTC(ii),par_ts.data(ii),tsamplUTC);
    end
    
    % Real mean 1-h and 3-h PAR if time series has finer measurement period (tm)
    % 3h follows 0:3:21 fixed grid. 1h mean is hour prior to sampling
    if tm <= 3
        tvecUTC = datevec(par_ts.mtimeUTC(iref_all));
        tsvecUTC = datevec(tsamplUTC);
        t3h = 0:3:21;
        tdif3h = t3h - tsvecUTC(4);
        t3h = t3h(tdif3h>-3 & tdif3h<=0);
        i3hSAT = iref_all(tvecUTC(:,4) >= t3h  & tvecUTC(:,4) < (t3h+3));
        if sum(i3hSAT)
            par0plus.tsam3hSAT = nanmean(par_ts.data(i3hSAT));
            par0plus.Ntsam3hSAT = sum(~isnan(i3hSAT));
        end
        i3h = par_ts.mtimeUTC(iref_all) > (tsamplUTC-3/24) & par_ts.mtimeUTC(iref_all) <= tsamplUTC;
        if sum(i3h)
            par0plus.tsam3h = nanmean(par_ts.data(iref_all(i3h)));
            par0plus.Ntsam3h = sum(i3h);
        end
        if tm <= 1
            i1h = par_ts.mtimeUTC(iref_all) > (tsamplUTC-1/24) & par_ts.mtimeUTC(iref_all) <= tsamplUTC;
            if sum(i1h)
                par0plus.tsam1h = nanmean(par_ts.data(iref_all(i1h)));
                par0plus.Ntsam1h = sum(i1h);
            end
        end
    end
    
    % 2.1 Mean PAR during 24 hours prior to sampling, corresponding data counts
    p24h = 24;
    i24h = (par_ts.mtimeUTC <= par_ts.mtimeUTC(iref)) & (par_ts.mtimeUTC > par_ts.mtimeUTC(iref) - p24h/24);
    if sum(i24h)
        par0plus.p24h = nanmean(par_ts.data(i24h));
        par0plus.Np24h = sum(i24h);
        par0plus.p24hmax = nanmax(par_ts.data(i24h));
    end
    
    % 2.2 Mean PAR during X hours prior to sampling, corresponding data counts
    if ~isnan(pXh) && ~isempty(pXh)
        iXh = (par_ts.mtimeUTC <= par_ts.mtimeUTC(iref)) & (par_ts.mtimeUTC > par_ts.mtimeUTC(iref) - pXh/24);
        iXh(isnan(iXh)) = [];
        if sum(iXh)
            par0plus.pXh = nanmean(par_ts.data(iXh));
            par0plus.NpXh = sum(iXh);
        end
    end
end

% LOCAL DAY DATA
% Note that local time is not exactly centered around solar noon because it
% follows time zone convention
tsvecUTC = datevec(tsamplUTC);
UTC_noon = solar_noon(tsvecUTC(1),yearday(tsvecUTC(3),tsvecUTC(2),tsvecUTC(1),0,0,0),24*(tsamplUTC-floor(tsamplUTC)),lon);
real_diff = UTC_noon - 12; % more exact but almost equivalent to [real_diff = -lon/15]
par_ts.mtimeLOC = par_ts.mtimeUTC - UTC_diff/24;
tvecLOC = datevec(par_ts.mtimeLOC);
tsamplLOC = tsamplUTC - UTC_diff/24;
tsvecLOC = datevec(tsamplLOC);
iday = (tvecLOC(:,1) == tsvecLOC(1)) & (tvecLOC(:,2) == tsvecLOC(2)) & (tvecLOC(:,3) == tsvecLOC(3));

if sum(iday)
    
    % 3. Mean daily PAR during the LOCAL day, data counts
    data_iday = par_ts.data(iday);
    par0plus.dayLOC = nanmean(data_iday);
    par0plus.NdayLOC = sum(iday);
    par0plus.dayLOCmax = nanmax(data_iday);
    
    % 4. Mean PAR at solar noon on sampling day
    % Implement interpolations based on interp_PAR_tests.m
    
    % 4.1 Take value from closest measurement period (backwards temporal search
    % since we look at integration that started before noon time, not after)
    mtimeUTC_iday = par_ts.mtimeUTC(iday);
    mtimeUTC_iday_vec = datevec(mtimeUTC_iday);
    hourdiff = mtimeUTC_iday_vec(:,4) + mtimeUTC_iday_vec(:,5)/60 - UTC_noon;
    mhd = nanmax(hourdiff(hourdiff<=0));
    if ~isnan(mhd) && (abs(mhd)/24 < tsearch)
        par0plus.noonLOC = data_iday(hourdiff==mhd);
    end
    % 3h follows 0:3:21 fixed grid. 1h mean is hour prior to sampling
    if tm <= 3
        mtimeLOC_iday = par_ts.mtimeLOC(iday);
        mtimeLOC_iday_vec = datevec(mtimeLOC_iday);
        if real_diff < UTC_diff
            inoon3h = mtimeLOC_iday_vec(:,4) >= 9 & mtimeLOC_iday_vec(:,4) < 12;
            if sum(inoon3h)
                par0plus.noon3hSATLOC = nanmean(data_iday(inoon3h));
                par0plus.Nnoon3hSATLOC = sum(inoon3h);
            end
        else
            inoon3h = mtimeLOC_iday_vec(:,4) >= 12 & mtimeLOC_iday_vec(:,4) < 15;
            if sum(inoon3h)
                par0plus.noon3hSATLOC = nanmean(data_iday(inoon3h));
                par0plus.Nnoon3hSATLOC = sum(inoon3h);
            end
        end
    end
    if tm <= 1
        inoon1h = hourdiff > -0.5 & hourdiff <= 0.5;
        if sum(inoon1h)
            par0plus.noon1hLOC = nanmean(data_iday(inoon1h));
            par0plus.Nnoon1hLOC = sum(inoon1h);
        end
    end
    
    % 4.2 Interpolation method 1: pchip --> Not implemented because it does not
    % provide good results fr tie series with less than 1 datum per hour, as
    % shown in code "interp_PAR_tests.m"
    
    % 4.3 Interpolation method 2: spline --> Not implemented because it does not
    % provide good results fr tie series with less than 1 datum per hour, as
    % shown in code "interp_PAR_tests.m"
    
end

% UTC DAY DATA
tvecUTC = datevec(par_ts.mtimeUTC);
idayUTC = (tvecUTC(:,1) == tsvecUTC(1)) & (tvecUTC(:,2) == tsvecUTC(2)) & (tvecUTC(:,3) == tsvecUTC(3));

if sum(idayUTC)
    % 3. Mean daily PAR during the UTC day, data counts
    data_idayUTC = par_ts.data(idayUTC);
    par0plus.dayUTC = nanmean(data_idayUTC);
    par0plus.NdayUTC = sum(idayUTC);
    % par0plus.dayLOCmaxUTC = nanmax(data_idayUTC);
end
