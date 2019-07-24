% PLOT DMS PROFILES FROM GREEN EDGE
clear
clc
% close all

load gcms_greenedge_proc
load ~/Desktop/GreenEdge/Irradiance/CASTS_ICE_matched_25-Apr-2017.mat
ICE = OUT; headICE = header_out; clear OUT header_out
load ~/Desktop/GreenEdge/meteo/GreenEdge_MET_2016_MGT.mat
load ~/Desktop/GreenEdge/meteo/indexs_METEO_c24h_p24h_p72h.mat
load ~/Desktop/GreenEdge/CTD_rosette/GreenEdge_CTD_2016_MGT.mat
path_anp = '~/Desktop/GreenEdge/Achim_Randelhoff/data/Randelhoff-et-al_GreenEdge_ANP_v0.1.csv';
ANP = csvread(path_anp, 1, 1);
headANP = {'station' 'depth' 'anp'};
ANP(ANP(:,strcmp(headANP,'depth'))>100,:) = [];

%% INITIAL SETTINGS
idms = 46;
idmspt = 49;
zSURF = 7;

%% PREPARE DATA

DATA = [data DMSOUT];
DATA_CTD = [];
headctd = (CTD{1}.head)';
head_data = {'year' 'month' 'day' 'nouse' 'nouse' 'stn' 'cast' 'nouse' 'nouse'...
    'nouse' 'nouse' 'niskin' 'depth' 'vol' 'prescreen'}';
headDATA_CTD = [head_data; repmat({'nouse'},size(data,2)-length(head_data),1);...
    headout; ({'temp' 'sal' 'sigt' 'N2' 'O2' 'CDOM' 'NO3' 'cp' 'cpsmooth1' 'anp'})'];
SURF = [];
headSURF = ({'stn' 'lat' 'lon' 'year' 'month' 'day' 'dms'...
    'SICm2d' 'SICm1d' 'SICday' 'wsp72' 'wsp24' 'wsc24'...
    'sst' 'sal' 'fdmsW97c24'  'fdmsN00c24' 'sigt' 'O2' 'CDOM' 'NO3' 'cpsmooth1'...
    'mld03' 'mld125' 'N2max03' 'zN2max03'  'N2max125' 'zN2max125' 'dmspt' 'dbm'})';
clear DMSOUT data

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
% xmax2 = [10 30 30 10];
xmax2 = [24 75 46 12];
zmax = [62 72 100 80];
lat = [68 70 70.5 69.5];
lon = {[-62.42 -61.6 60.82 -60 -58.93 -59.2 -57.77],...
    [-59.1 -58.66 -60.35 -61.25 -62.42],...
    [-64 -63.05 -62.42 -59.5 -60.42 -62.63 -61.62],...
    [-58.72 -57.87 -59.8 -61.58 -62.43 -63.23]};


%% CREATE SURFACE DATASET (USING LOOP OVER TRANSECTS)
% Match meteo time series, calculate DMS flux with minute resolution, then average over periods
% Add also MLD, sigma-t in the mixed layer, maximum N2 (B-V freq) between
% base of MLD and 100 m, and corresponding depth

for j = 1:length(lon) % loop on transects
    
    jlons = lon{j};
    jtrans = transects{j};
    jlat = lat(j);
    
    % Sort stations by longitude
    l_t = [jlons' jtrans']; l_t = sortrows(l_t,1);
    jlons = l_t(:,1); jtrans = l_t(:,2);
    
    for k = 1:length(jlons) % loop on transect stations
        
        kstn = jtrans(k);
        tmpall = DATA(station==kstn,:);
        anp.anp = ANP(ANP(:,strcmp(headANP,'station'))==kstn,strcmp(headANP,'anp'));
        anp.depth = ANP(ANP(:,strcmp(headANP,'station'))==kstn,strcmp(headANP,'depth'));
        % Remove if (DMS and DMSPt) or cast number absent
        tmpall((isnan(tmpall(:,idms)) & isnan(tmpall(:,idmspt))) | isnan(tmpall(:,7)),:) = [];
        
        if ~isempty(tmpall)
            
            % Match casts
            ncast = mode(tmpall(:,7));
            castmatch = ICE(:,strcmp(headICE,'cast')) == ncast;
            tmp2 = ICE(castmatch,12:14); % Ice
            tmp3 = CTD{castmatch};
            
            % Extract CTD variables
            depth = tmp3.DATA(:,strcmp(headctd,'Pres'));
            temp = tmp3.DATA(:,strcmp(headctd,'Temp'));
            sal = tmp3.DATA(:,strcmp(headctd,'Sal'));
            sigt = tmp3.DATA(:,strcmp(headctd,'sigt'));
            N2 = tmp3.DATA(:,strcmp(headctd,'N2'));
            Trans = tmp3.DATA(:,strcmp(headctd,'Trans'));
            O2 = tmp3.DATA(:,strcmp(headctd,'O2'));
            CDOM = tmp3.DATA(:,strcmp(headctd,'CDOM'));
            NO3 = tmp3.DATA(:,strcmp(headctd,'NO3'));
            
            % Calculate cp from Transmittance profile
            cp_out = get_cp(Trans,depth,kstn);
            cp = cp_out.q;
            cpsmooth1 = cp_out.s;
            dbm = cp_out.zdbm;
            
            % Fill SURF matrix
            tmp = tmpall(tmpall(:,13)>zSURF,:);
            
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
                ictd = nanmin(find(~isnan(depth) & ~isnan(temp) & ~isnan(sal)));
                if ~isempty(ictd) && depth(ictd)<=zSURF
                    APPEND(14:15) = [temp(ictd) sal(ictd)];
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
                
                % Append more CTD data
                if ~isempty(ictd) && depth(ictd)<=zSURF
                    APPEND(18:22) = [sigt(ictd) O2(ictd) CDOM(ictd) NO3(ictd) cpsmooth1(ictd)];
                end
                
                % MLD
                [mld,N2max,zN2max] = get_mld(depth,sigt,N2,kstn);
                APPEND(23:28) = [mld.d03 mld.d125 N2max.d03 N2max.d125 zN2max.d03 zN2max.d125];
                
                % DMSPt
                APPEND(1,29) = nanmean(tmp(:,idmspt),1);
                
                % Depth of cp maximum
                APPEND(1,30) = dbm;
                
                % Append SURF matrix
                SURF = [SURF; APPEND];
                
            end
            
            % Match and fill DATA_CTD matrix
            CTDMACTH = [temp sal sigt N2 O2 CDOM NO3 cp cpsmooth1];
            
            % Loop over profile depths and match CTD, ANP (which has
            % different depth axis)
            APPEND2 = [];
            append3 = [];
            for iz = 1:size(tmpall,1)
                zint = 0.5;
                zmeas = tmpall(iz,strcmp(headDATA_CTD,'depth'));
                if zmeas == 0.7, zint = 1.5; end % ensure match for surface data
                ictdmatch = depth<zmeas+zint & depth>zmeas-zint;
                if sum(ictdmatch)
                    APPEND2 = [APPEND2; nanmean(CTDMACTH(ictdmatch,:),1)];
                else
                    APPEND2 = [APPEND2; nan(1,size(CTDMACTH,2))];
                end
                ianpmatch = anp.depth<zmeas+zint & anp.depth>zmeas-zint;
                if sum(ictdmatch)
                    append3 = [append3; nanmean(anp.anp(ianpmatch,:),1)];
                else
                    append3 = [append3; nan];
                end
            end
            
            % Append DATA_CTD matrix
            DATA_CTD = [DATA_CTD; [tmpall APPEND2 append3]];
            
        end
    end
end

%% Save data

save('surf_profile_ALL_4plot.mat','SURF','headSURF','DATA_CTD','headDATA_CTD')

