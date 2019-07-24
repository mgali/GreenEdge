% INTEGRATE PAR FROM SBDART FOR ICE CAMP TIME SERIES

%% 2015

% Change to v3 as soon as available? finally v1 consensus is the best for
% ice camp
datadir = '~/Desktop/GreenEdge/Irradiance/SBDART_LUTs_outputs/Ed0plus_MODISA_LUT_SurfAlb_v1_1h_consensus';
year = 2016;
hperiod = 1;

% Spectral range and resolution
specres = 5;
wl = 290:specres:700;
ipar = wl >= 400 & wl <= 700;
iuvr = wl >= 290 & wl <= 400;

% Time axis
xmin = 119;
xmax = 189;
doys = xmin:xmax;
dd = (xmin+hperiod/2/24):hperiod/24:(xmax+1);
par_ts.mtimeUTC =  (datenum(2015,0,dd,0,0,0))';

% Coordinates
lat = 67.480;
lon = -63.790;

% Fill in loop
ED0 = nan(1,length(wl)); % length(doys)*(24/hperiod)
for j = 1:length(doys)
    TMP = dlmread(sprintf('%s/Ed0_%04i_%03i_%0.3f_%0.3f.txt',datadir,year,doys(j),lat,lon));
    ED0 = [ED0; TMP'];
end
ED0(1,:) = [];

%% Integrate spectrum
par_ts.data = mean(ED0(:,ipar),2)*(sum(ipar)-1)*specres;
par_ts.dataUVR = mean(ED0(:,iuvr),2)*(sum(iuvr)-1)*specres;

%% Save
save(sprintf('~/Desktop/GreenEdge/Irradiance/SBDART_LUTs_outputs/PAR_SBDART_GE%04i-ICECAMP_SurfAlb_v1_1h_consensus.mat',year))