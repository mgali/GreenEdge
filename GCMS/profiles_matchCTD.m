% MATCH DMS PROFILES TO CTD DATA
clear
clc
close all


%% INITIAL SETTINGS
% Rename data matrices and headers for clarity

% DMS data
load gcms_greenedge_proc
head_rawdms = {'year' 'month' 'day' 'nouse' 'nouse' 'stn' 'cast' 'nouse' 'nouse'...
    'nouse' 'nouse' 'niskin' 'depth' 'vol' 'prescreen'}';
head_rawdms = [head_rawdms; repmat({'nouse'},size(data,2)-length(head_rawdms),1)];
DATA = [data DMSOUT];
head_data = [head_rawdms; headout];
clear DMSOUT data headout head_rawdms

% CTD data
load ~/Desktop/GreenEdge/CTD_rosette/GreenEdge_CTD_2016_MGT.mat
idms = 46;
idmspt = 49;
head_ctd = (CTD{1}.head)';

%% Surface data (ice, clean meteo, meteo indexs for each cast)
load ~/Desktop/GreenEdge/Irradiance/CASTS_ICE_matched_25-Apr-2017.mat
ICE = OUT; head_ice = header_out; clear OUT header_out
load ~/Desktop/GreenEdge/meteo/GreenEdge_MET_2016_MGT.mat
load ~/Desktop/GreenEdge/meteo/indexs_METEO_c24h_p24h_p72h.mat


%% PREPARE DATA

% Create station list
station = DATA(:,6);
diffstation = diff(station);
stnax = [station(diffstation~=0); station(end)];
stnaxstr = cellstr(num2str(stnax));

% Define transects and longitudes
transects = {[400 403 406 409 413 412 418],...
    [507 506 512 515 519],...
    [600 603 605 615 612 604.5 608],...
    [703 700 707 713 716 719]};
lat = [68 70 70.5 69.5];
lon = {[-62.42 -61.6 60.82 -60 -58.93 -59.2 -57.77],...
    [-59.1 -58.66 -60.35 -61.25 -62.42],...
    [-64 -63.05 -62.42 -59.5 -60.42 -62.63 -61.62],...
    [-58.72 -57.87 -59.8 -61.58 -62.43 -63.23]};


%% MACTH CTD DATA

zSURF = 7;
SURF = [];
headSURF = {'stn' 'lat' 'lon' 'year' 'month' 'day' 'dms' 'SICm2d' 'SICm1d' 'SICday' 'wsp72' 'wsp24' 'wsc24' 'sst' 'sal' 'fdmsW97c24'  'fdmsN00c24'};

for j = 1:length(lon) % loop on transects
    
    jlons = lon{j};
    jtrans = transects{j};
    jlat = lat(j);
    
    % Sort stations by longitude
    l_t = [jlons' jtrans']; l_t = sortrows(l_t,1);
    jlons = l_t(:,1); jtrans = l_t(:,2);
    
    for k = 1:length(jlons) % loop on transect stations
        
        kstn = jtrans(k);
        tmp = DATA(station==kstn,:);
        % Remove if z not surface or DMS or cast number absent
        tmp(tmp(:,13)>zSURF | isnan(tmp(:,idms)) | isnan(tmp(:,7)),:) = [];
        
        % Match casts
        ncast = mode(tmp(:,7));
        castmatch = ICE(:,strcmp(head_ice,'cast')) == ncast;
        tmp2 = ICE(castmatch,12:14); % Ice
        tmp3 = CTD{castmatch};
        
        % Fill SURF matrix
        if ~isempty(tmp)
            APPEND = nan(1,length(headSURF));
            APPEND(1,1) = kstn; % station number
            APPEND(1,2:3) = [jlat jlons(k)]; % lat lon
            APPEND(1,4:7) = nanmean(tmp(:,[1:3 idms]),1); % date and DMS
            APPEND(1,8:10) = tmp2; % Ice
            
            % Macth meteo data
            imet.p72h = IMET.p72h{castmatch};
            imet.p24h = IMET.p24h{castmatch};
            imet.c24h = IMET.c24h{castmatch};
            APPEND(11:13) = [nanmean(met.ws(imet.p72h)) nanmean(met.ws(imet.p24h)) nanmean(met.ws(imet.c24h))];
            
            % Match CTD data
            depth = tmp3.DATA(:,strcmp(head_ctd,'Pres'));
            temp = tmp3.DATA(:,strcmp(head_ctd,'Temp'));
            sal = tmp3.DATA(:,strcmp(head_ctd,'Sal'));
            ictd = nanmin(find(~isnan(depth) & ~isnan(temp) & ~isnan(sal)));
            if ~isempty(ictd) && depth(ictd)<=zSURF
                APPEND(14:15) = [temp(ictd) sal(ictd)];
            else
                APPEND(14:15) = [nan nan];
            end
            
            % Calculate DMS flux over 24h centered on sampling (c24)
            finput.ws = met.ws(imet.c24h);
            finput.dms = APPEND(7)*ones(size(finput.ws));
            finput.sst = APPEND(14)*ones(size(finput.ws));
            finput.sal = APPEND(15)*ones(size(finput.ws));
            if isnan(APPEND(14)), finput.sst = -0.2402; end % assign mean SURF sst if not available
            if isnan(APPEND(15)), finput.sal = 32.0936; end % assign mean SURF sal if not available
            fout.W97 = fdms(finput.dms,finput.ws,finput.sst,finput.sal,'W97');
            fout.N00 = fdms(finput.dms,finput.ws,finput.sst,finput.sal,'N00GBC');
            APPEND(16:17) = [nanmean(fout.W97) nanmean(fout.N00)];
            
            % Append SURF matrix
            SURF = [SURF; APPEND];
            
        end
        
        
    end
end



%% Classify by ice cover (open receding covered)

icemean = nanmean(SURF(:,8:10),2);
icemin = nanmin(SURF(:,8:10),[],2);
icemax = nanmax(SURF(:,8:10),[],2);
icecrit1 = 0.10;
icecrit2 = 0.80;
icebin = 2*ones(size(icemax));
icebin(icemax<icecrit1) = 1; % Open: Persistent <10% ice
icebin(icemin>icecrit2) = 3; % Covered: Persistent >80% ice
corrflux = SURF(:,16).*(1-icemean);

%% Boxplot of DMS by ice cover

% xtl = {sprintf('Open (N=%i)',sum(icebin==1));
%     sprintf('Receding (N=%i)',sum(icebin==2));
%     sprintf('Ice (N=%i)',sum(icebin==3))};
xtl = {sprintf('Open: SIC<%0.2f',icecrit1);
    sprintf('Receding');
    sprintf('Ice: SIC>%0.2f',icecrit2)};
ymin = -1;
ymax = 27.5;
cc = {'r' 'k' 'b'};

figure(1), clf

subplot(1,2,1)
boxplot(SURF(:,7),icebin,'colors','rkb','medianstyle','line'), hold on
for ii = 1:3
    plot(ii,nanmean(SURF(icebin==ii,7)),'.','markersize',30,'color',cc{ii})
    text(ii-0.1,ymax+0.5,sprintf('N=%i',sum(icebin==ii)),'fontsize',18)
end
ylabel('DMS (nM)','fontsize',18)
set(gca,'tickdir','out','xtick',1:max(icebin),'xticklabel',xtl,'fontsize',18)
ylim([ymin ymax])

subplot(1,2,2)
boxplot(corrflux,icebin,'colors','rkb','medianstyle','line'), hold on
cc = {'r' 'k' 'b'};
for ii = 1:3
    plot(ii,nanmean(corrflux(icebin==ii)),'.','markersize',30,'color',cc{ii})
    text(ii-0.1,ymax+0.5,sprintf('N=%i',sum(icebin==ii)),'fontsize',18)
end
ylabel('F_{DMS} (µmol m^{-2} d^{-1})','fontsize',18)
set(gca,'tickdir','out','xtick',1:max(icebin),'xticklabel',xtl,'fontsize',18)
ylim([ymin ymax])
