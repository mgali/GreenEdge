% ALLSURF SURFACE DMS AND ADDITIONAL CTD-DERIVED DATA TO ACHIM'S STATION
% MATRIX DATA
clear
clc

load samples.GE2016only.castsOnly.mat
GE2016 = samplesGE2016castsOnly;
load gcms_greenedge_proc
load ~/Desktop/GreenEdge/CTD_rosette/GreenEdge_CTD_2016_MGT.mat
headCTD = CTD{1}.head;
load surf_profile_ALL_4plot
path_rand = '~/Desktop/GreenEdge/Achim_Randelhoff/data/Randelhoff-et-al_Ice-edge-paper_per-station_v0.1.csv';
RAND = csvread(path_rand, 1, 0);
headRAND = ({'station' 'isolume_m_at_0.1_Einm-2d-1' 'isolume_m_at_0.415_Einm2-d1'...
    'IceConcentration_percent' 'MixedLayerDepth_0.1gkg-1-criterion_m' 'hBD_m'...
    'Nitracline_m' 'AW-ArW_clustering_coefficient' 'OWD' 'PAR_at_3m_Einm-2d-1'})';
rst = RAND(:,strcmp('station',headRAND));

%% Load HPLC
headHPLC = {'date_anal' 'cruise' 'leg' 'type' 'lat' 'lon' 'date_sam' 't_sam' 'stn' 'code_sam' 'cast' 'niskin' 'depth' 'filt_vol(L)'...
    'chlc3' 'Chlc3-QA' 'chlc2group' 'Chlc2-QA' 'chldaSUM' 'Chlda-QA' 'peri' 'Peri-QA' 'phdaSUM'	'Phda-QA'...
    'uriol' 'Uriol-QA'	'but' 'But-QA' 'fuco' 'Fuco-QA' 'neo' 'Neo-QA' 'pras' 'Pras-QA'...
    'viola' 'viola-QA' 'hex' 'Hex-QA' 'asta' 'Asta-QA' 'micral' 'micral-QA' 'diadino' 'Diadino-QA'...
    'anthera' 'Anthera-QA' 'allo' 'Allo-QA' 'diato' 'Diato-QA' 'zea' 'Zea-QA' 'lut' 'Lut-QA' 'bchla' 'Bchla-QA'...
    'chlb' 'Chlb-QA' 'dvchla' 'DVChla-QA' 'chla' 'Chla-QA' 'tchla' 'Tchla-QA' 'phytnaSUM' 'Phytna-QA'...
    'caro_like_Prasi' 'caro_like_Prasi-QA' 'tcar' 'Tcar-QA' '19BF_like' '19BF_like-QA' '19HF_likeSUM' '19HF_likeSUM-QA'}';
HPLC = xlsread('~/Desktop/GreenEdge/GreenEdge-Amundsen-pigments-180131.xlsx');
z.HPLC = HPLC(:,strcmp('depth',headHPLC));
z.HPLC(z.HPLC==0) = 0.7;
HPLC(:,strcmp('depth',headHPLC)) = z.HPLC;
ni.HPLC = HPLC(:,strcmp('niskin',headHPLC));
ca.HPLC = HPLC(:,strcmp('cast',headHPLC));
ca.HPLC = ca.HPLC - 1601000;
HPLC(:,strcmp('cast',headHPLC)) = ca.HPLC;
chl.HPLC = HPLC(:,strcmp(headHPLC,'tchla'));

%% Check match
% figure(),plot(86:203,GE2016(86:203,end),'ok'), hold on
% plot(DATA_CTD(:,7),DATA_CTD(:,6),'.r','markersize',6)
% set(gca,'xtick',80:10:200)
% grid on

%% Apend data in loop

headALLSURF = [headRAND;...
    headSURF([8:15 18:22 23:24 25 27 26 28 30 29 7 16 17]);...
    {'idms_mld'; 'idms_hBD'; 'idms_isolu'; 'idms_z60';'idms_z10';...
    'idmspt_mld'; 'idmspt_hBD'; 'idmspt_isolu'; 'idmspt_z60'; 'idmspt_z10';...
    'icp_mld'; 'icp_hBD'; 'icp_isolu'; 'icp_z60'; 'icp_z10';...
    'Tchla'; 'iTchla_mld'; 'iTchla_hBD'; 'iTchla_isolu'; 'iTchla_z60'; 'iTchla_z10';...
    }];
to_correct = headALLSURF(26:29);
headALLSURF([27 29 26 28]) = to_correct;
ALLSURF = nan(size(GE2016,1),length(headALLSURF));

profstart = 86;
profend = size(GE2016,1);
for k = profstart:profend % loop on transect stations
    
    kstn = GE2016(:,strcmp(headerGE2016,'stn'));
    kk = kstn(k);
    kcast = k;
    fprintf('Station = %i, cast = %i\n',kk,k)
    
    % Fill in RAND data
    tmp1 = RAND(RAND(:,strcmp(headRAND,'station'))==kk,:);
    if ~isempty(tmp1)
        ALLSURF(k, 1:length(headRAND)) = tmp1;
    end
    
    % Fill CTD data irrespective of DMS being available
    
    % Extract CTD variables
    tmp3 = CTD{kcast};
    if isstruct(tmp3) && numel(fieldnames(tmp3))==2
        depth = tmp3.DATA(:,strcmp(headCTD,'Pres'));
        temp = tmp3.DATA(:,strcmp(headCTD,'Temp'));
        sal = tmp3.DATA(:,strcmp(headCTD,'Sal'));
        sigt = tmp3.DATA(:,strcmp(headCTD,'sigt'));
        N2 = tmp3.DATA(:,strcmp(headCTD,'N2'));
        Trans = tmp3.DATA(:,strcmp(headCTD,'Trans'));
        O2 = tmp3.DATA(:,strcmp(headCTD,'O2'));
        CDOM = tmp3.DATA(:,strcmp(headCTD,'CDOM'));
        NO3 = tmp3.DATA(:,strcmp(headCTD,'NO3'));
        
        % Process CTD profile and append mld and stratification data
        [mld,N2max,zN2max] = get_mld(depth,sigt,N2,kstn);
        ALLSURF(k,24:29) = [mld.d03 mld.d125 N2max.d03 N2max.d125 zN2max.d03 zN2max.d125];
        
        % Calculate cp profile and append DBM depth
        cp_out = get_cp(Trans,depth,kstn);
        cp = cp_out.q;
        ALLSURF(k,30) = cp_out.zdbm;
        
        % Add surface CTD data
        isf = depth<=2;
        ALLSURF(k,17:22) = nanmean([temp(isf) sal(isf) sigt(isf) O2(isf) CDOM(isf) NO3(isf)],1);
    end
    
    % Define bounds for vertical integrals, integrate
    zbounds.mld = round([0 ALLSURF(k,strcmp(headALLSURF,'MixedLayerDepth_0.1gkg-1-criterion_m'))]);
    zbounds.hBD = round([0 ALLSURF(k,strcmp(headALLSURF,'hBD_m'))]);
    zbounds.isolu = round([0 ALLSURF(k,strcmp(headALLSURF,'isolume_m_at_0.415_Einm2-d1'))]);
    zbounds.z60 = [0 60];
    zbounds.z10 = [0 10];
    
    if ~isnan(kk) && kk >= 400
        % Append surface and sea-air flux data
        tmpsurf = SURF(SURF(:,strcmp(headSURF,'stn'))==kk,:);
        if ~isempty(tmpsurf)
            ALLSURF(k,[11:16 31:34]) = tmpsurf([8:13 29 7 16 17]);
        else
            ALLSURF(k,[11:16 31:34]) = nan(1,10);
        end
        
        % DMS(P) profile data: vertical integrals
        tmpall = DATA_CTD(DATA_CTD(:,strcmp(headDATA_CTD,'stn'))==kk,:);
        if ~isempty(tmpall)
            z_in = tmpall(:,strcmp(headDATA_CTD,'depth'));
            % DMS
            mycol = strcmp(headDATA_CTD,'dms_consens_cf68');
            ALLSURF(k,35) = get_zintegral(z_in, tmpall(:,mycol), zbounds.mld, kk, kcast, 'dms', 'mld');
            ALLSURF(k,36) = get_zintegral(z_in, tmpall(:,mycol), zbounds.hBD, kk, kcast, 'dms', 'hBD');
            ALLSURF(k,37) = get_zintegral(z_in, tmpall(:,mycol), zbounds.isolu, kk, kcast, 'dms', 'isolu');
            ALLSURF(k,38) = get_zintegral(z_in, tmpall(:,mycol), zbounds.z60, kk, kcast, 'dms', 'z60');
            ALLSURF(k,39) = get_zintegral(z_in, tmpall(:,mycol), zbounds.z10, kk, kcast, 'dms', 'z10');
            % DMSPt
            mycol = strcmp(headDATA_CTD,'dmspt');
            ALLSURF(k,40) = get_zintegral(z_in, tmpall(:,mycol), zbounds.mld, kk, kcast, 'dmspt', 'mld');
            ALLSURF(k,41) = get_zintegral(z_in, tmpall(:,mycol), zbounds.hBD, kk, kcast, 'dmspt', 'hBD');
            ALLSURF(k,42) = get_zintegral(z_in, tmpall(:,mycol), zbounds.isolu, kk, kcast, 'dmspt', 'isolu');
            ALLSURF(k,43) = get_zintegral(z_in, tmpall(:,mycol), zbounds.z60, kk, kcast, 'dmspt', 'z60');
            ALLSURF(k,44) = get_zintegral(z_in, tmpall(:,mycol), zbounds.z10, kk, kcast, 'dmspt', 'z10');
        end % test non-empty subset
        
        % TChla (HPLC)
        tmp2 = HPLC(ca.HPLC==kcast,strcmp(headHPLC,'depth')|strcmp(headHPLC,'tchla'));
        if ~isempty(tmp2)
            tchla = tmp2(:,2);
            zchla = tmp2(:,1);
            ALLSURF(k,51) = nanmean(tchla(zchla>=2),1);
            ALLSURF(k,51) = get_zintegral(zchla, tchla, zbounds.mld, kk, kcast, 'tchla', 'mld');
            ALLSURF(k,52) = get_zintegral(zchla, tchla, zbounds.hBD, kk, kcast, 'tchla', 'hBD');
            ALLSURF(k,53) = get_zintegral(zchla, tchla, zbounds.isolu, kk, kcast, 'tchla', 'isolu');
            ALLSURF(k,54) = get_zintegral(zchla, tchla, zbounds.z60, kk, kcast, 'tchla', 'z60');
            ALLSURF(k,55) = get_zintegral(zchla, tchla, zbounds.z10, kk, kcast, 'tchla', 'z10');
        end % test non-empty subset
        
    end % test non-nan station number
    
    % cp
    if isstruct(tmp3) && numel(fieldnames(tmp3))==2
        ALLSURF(k,45) = get_zintegral(depth, cp_out.s, zbounds.mld, kk, kcast, 'cp', 'mld');
        ALLSURF(k,46) = get_zintegral(depth, cp_out.s, zbounds.hBD, kk, kcast, 'cp', 'hBD');
        ALLSURF(k,47) = get_zintegral(depth, cp_out.s, zbounds.isolu, kk, kcast, 'cp', 'isolu');
        ALLSURF(k,48) = get_zintegral(depth, cp_out.s, zbounds.z60, kk, kcast, 'cp', 'z60');
        ALLSURF(k,49) = get_zintegral(depth, cp_out.s, zbounds.z10, kk, kcast, 'cp', 'z10');
    end
    
end
save(sprintf('ALLSURF_%i_to_%i.mat',profstart,profend),'ALLSURF','headALLSURF')

%% Merge ALLSURF subsets
load ALLSURF_1_to_85; A1 = ALLSURF(1:85,:);
load ALLSURF_86_to_203; A2 = ALLSURF(86:203,:);
ALLSURF = [A1; A2];

%% Save ful ALLSURF
save('ALLSURF.mat','ALLSURF','headALLSURF')

