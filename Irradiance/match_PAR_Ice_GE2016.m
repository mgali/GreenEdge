% MATCH IRRADIANCE AND ICE DATA

%% General definitions
modays = [31 28 31 30 31 30 31 31 30 31 30 31];
modaysleap = [31 29 31 30 31 30 31 31 30 31 30 31];
mocum = cumsum(modays); mocum = [0 mocum(1:(end-1))];
mocumleap = cumsum(modaysleap); mocumleap = [0 mocumleap(1:(end-1))];

%% Load ship data

% ! head -1 station_data/LogBookScience_GreenEdge_leg1a.csv.txt
shipheader = ({'day' 'month' 'year' 'H' 'M' 'yearUTC' 'monthUTC' 'dayUTC' 'HUTC' 'MUTC' 'latNdegr'...
    'latNmin' 'lonWdegr' 'lonWmin' 'cap' 'cast' 'zbot' 'Wdir' 'Wspeed' 'Tair' 'Twater' 'Patm' 'RH' 'Ice'})';
leg1a = dlmread('~/Desktop/GreenEdge/Irradiance/station_data/LogBookScience_GreenEdge_leg1a.csv.nohead.txt');
leg1b = dlmread('~/Desktop/GreenEdge/Irradiance/station_data/LogBookScience_GreenEdge_leg1b.csv.nohead.txt');
leg1b(:,end) = [];
G = [leg1a; leg1b];

% Create cast list
lat = G(:,11) + G(:,12)/60;
lon = -(G(:,13) + G(:,14)/60);
month = G(:,7);
year = G(:,6);
day = G(:,8);
doy = nan(size(year,1),1);
isleap = ~mod(year,4);
doy(isleap) = (mocumleap(month(isleap)))' + day(isleap);
isnotleap = ~~mod(year,4);
doy(isnotleap) = (mocum(month(isnotleap)))' + day(isnotleap);
UTC_decimal = (G(:,9) + G(:,10)/60)/24;
cast = G(:,16);
zbot = G(:,17);
iceVisual = G(:,24);
dayc = datenum(year,month,day);

casts = [year month day doy UTC_decimal lat lon cast zbot iceVisual];
castsheader = ({'yearUTC' 'monthUTC' 'dayUTC' 'DOY' 'HUTC' 'lat' 'lon' 'cast' 'zbot' 'iceVisual'})';

%% Load ice data for day of sampling and previous 2 days
icefile = xlsread('~/Desktop/GreenEdge/Ice_Snow/GE2016/ice_history.xlsx','C2:I19857');

truedate1 = datenum(2016,6,9,0,0,0);
timelag1 = (icefile(1,1)) - truedate1;
truedate6 = datenum(2016,3,1,0,0,0);
timelag6 = (icefile(1,6)) - truedate6;
if timelag1==timelag6
    icefile(:,1) = icefile(:,1) - timelag1;
    icefile(:,6) = icefile(:,6) - timelag6;
else
    error('check dates again!')
end

% Reshape ice data
nicestn = 146;
ndays = size(icefile,1)/nicestn;
if ~mod(ndays,1)
%     LATI = reshape(icefile(:,2),nicestn,ndays);
%     LONI = reshape(icefile(:,3),nicestn,ndays);
%     DATESTNI = reshape(icefile(:,1),nicestn,ndays);
%     DATEI = reshape(icefile(:,end-1),nicestn,ndays);
    ICE = reshape(icefile(:,end),nicestn,ndays);
else
    error('check file!')
end

% % Make sure the list is always the same
% figure(),plot(DATEI)
% figure(),plot(DATEI')
% figure(),plot(DATESTNI')
% figure(),plot(DATESTNI')
% figure(),plot(DATESTNI)
% figure(),plot(LATI')
% figure(),plot(LATI)
% figure(),plot(LONI)

datestni = icefile(1:nicestn,1);
lati = icefile(1:nicestn,2);
loni = -icefile(1:nicestn,3);

% Make list of ice dates
kk = icefile(:,end-1);
datei = (kk(1):kk(end))';

%% Match by day and distance

ICEC = nan(size(casts,1),size(ICE,2));
dkm = nan(size(casts,1),1);

% Prealloacte ice concentration during day and n days backwards
ndhist = 4;
ICEC_HIST = nan(size(casts,1),ndhist+1);

castindexs = nan(size(ICE,1),1);

for ic = 1:size(casts,1)
    iday = dayc(ic) == datestni;
    if sum(iday)
        TMPICE = ICE(iday,:);
        [arclen,az] = distance(lat(ic),lon(ic),lati(iday),loni(iday));
        tmpdkm = arclen*60*1.852; % 60 nautical miles in one arcdegree, converted to km
        dkm(ic) = nanmin(tmpdkm);
        idkm = tmpdkm == dkm(ic);
        if sum(idkm)
            ICEC(ic,:) = TMPICE(idkm,:);
            
            % Find ice concentration during sampling day and n days backwards
            ihist = find(dayc(ic) == datei);
            ICEC_HIST(ic,:) = ICEC(ic,(ihist-ndhist):ihist);
            
        end
    end
end

OUT = [casts dkm ICEC_HIST];
icecheader = ({'SIC_minus2d' 'SIC_minus1d' 'SIC_day'})';
header_out = [castsheader; 'dkm'; icecheader];
save(sprintf('CASTS_ICE_matched_%s.mat',date),'OUT','header_out')

%% Compare visual ice concentration and satellite daily pixel data

% x = iceVisual;
% % x(x==0.5) = nan;
% y = ICEC_HIST(:,end);
% figure(4), clf
% plot(x(x~=0.5),y(x~=0.5),'ob','linewidth',2), hold on
% plot(x(x==0.5),y(x==0.5),'or','linewidth',2)
% xlabel('Sea-ice concentration during cast (visual)','fontsize',16)
% ylabel('Daily sea-ice concentration AMSR-2','fontsize',16)
% legend('All data','Berg (set to 0.5)','location','northoutside')
% plot([0 10],[0 1],'-k')
% [r,p] = corr(x(x~=0.5),y(x~=0.5),'rows','pairwise','type','pearson');
% [rs,ps] = corr(x(x~=0.5),y(x~=0.5),'rows','pairwise','type','spearman');
% text(0.5,0.9,sprintf('r = %0.2f, p = %0.0e',r,p),'color','b','fontsize',14)
% text(0.5,0.8,sprintf('r_{S} = %0.2f, p_{S} = %0.0e',rs,ps),'color','b','fontsize',14)
% set(gca,'fontsize',14)
% axis square
% box off

