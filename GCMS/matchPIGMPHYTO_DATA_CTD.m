clear, clc

%% Load data
load('surf_profile_ALL_4plot.mat','SURF','headSURF','DATA_CTD','headDATA_CTD')
clear headSURF; clear SURF;
compact = [1 2 3 6 7 12 13 45 48:54 56];
headDATA_CTD = headDATA_CTD(compact);
headDATA_CTD{8} = 'dms';
DATA_CTD = DATA_CTD(:,compact);
z.DATA_CTD = DATA_CTD(:,strcmp('depth',headDATA_CTD));
ni.DATA_CTD = DATA_CTD(:,strcmp('niskin',headDATA_CTD));
ca.DATA_CTD = DATA_CTD(:,strcmp('cast',headDATA_CTD));
st.DATA_CTD = DATA_CTD(:,strcmp('stn',headDATA_CTD));

%% Load turner pigments
headPITU = {'sampleID' 'date' 'stn' 'cast' 'WaterOrIce' 'niskin' 'depth' 'chla(ug/L)' 'phaeo(ug/L)'}';
iPITU = [8 9];
PITU = xlsread('~/Desktop/GreenEdge/GE_Amundsen_2016_ExtractedChl.xlsx');
z.PITU = PITU(:,strcmp('depth',headPITU));
z.PITU(z.PITU==1) = 0.7;
PITU(:,strcmp('depth',headPITU)) = z.PITU;
ni.PITU = PITU(:,strcmp('niskin',headPITU));
ca.PITU = PITU(:,strcmp('cast',headPITU));
st.PITU = PITU(:,strcmp('stn',headPITU));

%% Load HPLC
headHPLC = {'date_anal' 'cruise' 'leg' 'type' 'lat' 'lon' 'date_sam' 't_sam' 'stn' 'code_sam' 'cast' 'niskin' 'depth' 'filt_vol(L)'...
    'chlc3' 'Chlc3-QA' 'chlc2group' 'Chlc2-QA' 'chldaSUM' 'Chlda-QA' 'peri' 'Peri-QA' 'phdaSUM'	'Phda-QA'...
    'uriol' 'Uriol-QA'	'but' 'But-QA' 'fuco' 'Fuco-QA' 'neo' 'Neo-QA' 'pras' 'Pras-QA'...
    'viola' 'viola-QA' 'hex' 'Hex-QA' 'asta' 'Asta-QA' 'micral' 'micral-QA' 'diadino' 'Diadino-QA'...
    'anthera' 'Anthera-QA' 'allo' 'Allo-QA' 'diato' 'Diato-QA' 'zea' 'Zea-QA' 'lut' 'Lut-QA' 'bchla' 'Bchla-QA'...
    'chlb' 'Chlb-QA' 'dvchla' 'DVChla-QA' 'chla' 'Chla-QA' 'tchla' 'Tchla-QA' 'phytnaSUM' 'Phytna-QA'...
    'caro_like_Prasi' 'caro_like_Prasi-QA' 'tcar' 'Tcar-QA' '19BF_like' '19BF_like-QA' '19HF_likeSUM' '19HF_likeSUM-QA'}';
iHPLC = 15:74;
HPLC = xlsread('~/Desktop/GreenEdge/GreenEdge-Amundsen-pigments-180131.xlsx');
z.HPLC = HPLC(:,strcmp('depth',headHPLC));
z.HPLC(z.HPLC==0) = 0.7;
HPLC(:,strcmp('depth',headHPLC)) = z.HPLC;
ni.HPLC = HPLC(:,strcmp('niskin',headHPLC));
ca.HPLC = HPLC(:,strcmp('cast',headHPLC));
ca.HPLC = ca.HPLC - 1601000;
HPLC(:,strcmp('cast',headHPLC)) = ca.HPLC;
st.HPLC = HPLC(:,strcmp('stn',headHPLC));


%% MICROSCOPY TAXONOMY
[TAXO,headTAXO] = f_readxls_taxo('~/Desktop/GreenEdge','GE2016_AMD_Taxonomy-Microscopy_29082017_MG.xlsx');
iTAXO = 3:length(headTAXO);
z.TAXO = TAXO(:,strcmp('depth',headTAXO));
ca.TAXO = TAXO(:,strcmp('cast',headTAXO));

%% CYTOMETRY

%% IFCB

%% MATCH IN LOOP

% Preallocate
A.PITU = nan(size(DATA_CTD,1),length(iPITU));
A.HPLC = nan(size(DATA_CTD,1),length(iHPLC));
A.TAXO = nan(size(DATA_CTD,1),length(iTAXO));

dz = [0.5 1]; % Search radii in case no data found
dzflag.PITU = nan(size(A.PITU,1),1);
dzflag.HPLC = nan(size(A.HPLC,1),1);
dzflag.TAXO = nan(size(A.TAXO,1),1);

for j = 1:size(DATA_CTD,1)
    
    % Match Turner pigments using cast number and depth (niskin# is incomplete, do not use)
    mPITU = (z.PITU==z.DATA_CTD(j)) & (ca.PITU==ca.DATA_CTD(j));
    if sum(mPITU)
        dzflag.PITU(j) = 0; % sprintf('Matchup found for row %i!',j)
    else
        sprintf('Matchup not found for row %i!',j)
        mPITU = (z.PITU>=z.DATA_CTD(j)-dz(1) & z.PITU<=z.DATA_CTD(j)+dz(1)) & (ca.PITU==ca.DATA_CTD(j));
        if sum(mPITU)
            dzflag.PITU(j) = 1;
        else
            mPITU = (z.PITU>=z.DATA_CTD(j)-dz(2) & z.PITU<=z.DATA_CTD(j)+dz(2)) & (ca.PITU==ca.DATA_CTD(j)); 
            if sum(mPITU), dzflag.PITU(j) = 2; end
        end
    end
    A.PITU(j,:) = nanmean(PITU(mPITU,iPITU),1);
    
    % Match HPLC pigments using cast number and depth
    mHPLC = (z.HPLC==z.DATA_CTD(j)) & (ca.HPLC==ca.DATA_CTD(j));
    if sum(mHPLC)
        dzflag.HPLC(j) = 0; % sprintf('Matchup found for row %i!',j)
    else
        sprintf('Matchup not found for row %i!',j)
        mHPLC = (z.HPLC>=z.DATA_CTD(j)-dz(1) & z.HPLC<=z.DATA_CTD(j)+dz(1)) & (ca.HPLC==ca.DATA_CTD(j));
        if sum(mHPLC)
            dzflag.HPLC(j) = 1;
        else
            mHPLC = (z.HPLC>=z.DATA_CTD(j)-dz(2) & z.HPLC<=z.DATA_CTD(j)+dz(2)) & (ca.HPLC==ca.DATA_CTD(j)); 
            if sum(mHPLC), dzflag.HPLC(j) = 2; end
        end
    end
    A.HPLC(j,:) = nanmean(HPLC(mHPLC,iHPLC),1);
    
    % Match TAXO counts using cast number and depth
    mTAXO = (z.TAXO==z.DATA_CTD(j)) & (ca.TAXO==ca.DATA_CTD(j));
    if sum(mTAXO)
        dzflag.TAXO(j) = 0; % sprintf('Matchup found for row %i!',j)
    else
        sprintf('Matchup not found for row %i!',j)
        mTAXO = (z.TAXO>=z.DATA_CTD(j)-dz(1) & z.TAXO<=z.DATA_CTD(j)+dz(1)) & (ca.TAXO==ca.DATA_CTD(j));
        if sum(mTAXO)
            dzflag.TAXO(j) = 1;
        else
            mTAXO = (z.TAXO>=z.DATA_CTD(j)-dz(2) & z.TAXO<=z.DATA_CTD(j)+dz(2)) & (ca.TAXO==ca.DATA_CTD(j)); 
            if sum(mTAXO), dzflag.TAXO(j) = 2; end
        end
    end
    A.TAXO(j,:) = nanmean(TAXO(mTAXO,iTAXO),1);
    
end


%% APPEND AND SAVE
DATA = [DATA_CTD A.PITU A.HPLC];
headDATA = [headDATA_CTD; headPITU(iPITU); headHPLC(iHPLC)];


%% Quick scatterplot DMS vs. total carotenoids
% cmap = flip(parula);
% figure(11),
% ix = strcmp('tcar',headDATA);
% iy = strcmp('dms',headDATA);
% iz = strcmp('depth',headDATA);
% hold on
% scatter(DATA(:,ix),DATA(:,iy),50,DATA(:,iz),'fill'); colormap(cmap);
% title(sprintf('r = %0.2f',corr(DATA(:,ix),DATA(:,iy),'rows','pairwise')),'fontsize',16)
% xlabel(headDATA{ix},'fontsize',16)
% ylabel(headDATA{iy},'fontsize',16)
% c = colorbar;
% ylabel(c,headDATA{iz},'fontsize',16)


%% Quick correlation table and multilinear regression
ic = 19:2:size(DATA,2);
XC = DATA(:,ic);
removeX = sum(isnan(XC))>30;
XC(:,removeX) = [];
ic(removeX) = [];
yc = DATA(:,8);
R = corr(XC,yc,'rows','pairwise');
figure(2),
bar(R)
xlim([0 length(ic)+1])
set(gca,'tickdir','out','ticklength',[.005 .005],'xtick',1:length(ic),'xticklabel',[])
for j = 1:length(ic)
    text(j,-0.12,cellstr(headDATA(ic(j))),'rotation',-45,'interpreter','none')
end
title('Correlation coefficient: DMS vs. Pigments','fontsize',16)


%%
% removeX = sum(isnan(XC))>30;
% XC(:,removeX) = [];
% ic(removeX) = [];
remove = isnan(yc) | sum(isnan(XC),2)>0;
sum(~remove)
[b,bi,r,ri,st] = regress(yc(~remove),[ones(size(yc(~remove))) XC(~remove,:)]);
figure(2),
errorbar(1:length(ic)+1,b,bi(:,1),bi(:,2),'ok')
xlim([0 length(ic)+2])
set(gca,'tickdir','out','ticklength',[.005 .005],'xtick',1:length(ic),'xticklabel',[])
for j = 1:length(ic)+1
    if j == 1
        text(j,-0.22,'Intercept','rotation',-45,'interpreter','none')
    else
        text(j,-0.22,cellstr(headDATA(ic(j-1))),'rotation',-45,'interpreter','none')
    end
end
title(sprintf('Linear regresion DMS vs. Pigments, R^2 = %0.2f, p = %0.2e',st(1),st(3)),...
    'fontsize',16)


%% Quick photoprotection indices: deepoxidation and similar
% figure(12), clf
% ix1 = strcmp('chla',headDATA);
% ix2 = strcmp('diadino',headDATA);
% ix3 = strcmp('diato',headDATA);
% iy = strcmp('dms',headDATA);
% iz = strcmp('depth',headDATA);
% x = DATA(:,ix3)./(DATA(:,ix2)+DATA(:,ix3)); xlab = 'Dt/(Dd+Dt)';
% % x = (DATA(:,ix2)+DATA(:,ix3))./DATA(:,ix1); xlab = '(Dd+Dt)/Chla';
% y = DATA(:,iy); y(y>35)=nan;
% z = DATA(:,iz);
% hold on
% scatter(x,y,50,z,'fill'); colormap(cmap);
% title(sprintf('r = %0.2f',corr(x,y,'rows','pairwise')),'fontsize',16)
% xlabel(xlab,'fontsize',16)
% ylabel(headDATA{iy},'fontsize',16)
% c = colorbar;
% ylabel(c,headDATA{iz},'fontsize',16)


%% Quick correlation between DMS and TAXO
ic = iTAXO-2;
XC = A.TAXO(:,ic);
removeX = sum(XC==0)>5;
XC(:,removeX) = [];
ic(removeX) = [];
yc = DATA(:,8);
remove = DATA(:,strcmp(headDATA,'depth'))<=6; %TO SEPARATE SURFACE AND SCM
XC(remove,:) = nan;
yc(remove) = nan;
Rp = corr(XC,yc,'rows','pairwise');
Rs = corr(XC,yc,'rows','pairwise','type','Spearman');
figure(13),
bar(Rs,'barwidth',0.5,'facecolor',[.6 .6 .6]) % [Rp Rs]
xlim([0 length(ic)+1])
ylim([-0.7 1])
set(gca,'tickdir','out','ticklength',[.005 .005],'xtick',1:length(ic),...
'xticklabel',[],'fontsize',16,'position',[.15 .25 .7 .6])
ylim([-1 1])
for j = 1:length(ic)
    text(j,-1.05,cellstr(headTAXO(ic(j)+2)),'rotation',-45,'interpreter','none','fontsize',16)
end
ylabel('r (Spearman)','fontsize',16)
title('[DMS] vs. phyto. counts','fontsize',16)
% legend('r_P','r_S')


%% Quick boxplot of cell counts in major groups plus Phaeocystis
figure(14), clf
bar(R)
boxplot(log10(XC+1));
xlim([0 length(ic)+1])
ylim([1 7.2]);
set(gca,'tickdir','out','ticklength',[.005 .005],'xtick',[],'xticklabel',[])
for j = 1:length(ic)
    text(j,0.7,cellstr(headTAXO(ic(j)+2)),'rotation',-45,'interpreter','none')
end
title('Cells/L (microscopy)','fontsize',16)


%% Quick scatterplot between DMS and TAXO, normalized by...
cmap = flip(parula);
figure(15),
counts = nansum(A.TAXO(:,1:11),2); % total
diatom = nansum(A.TAXO(:,1:2),2); % diatom
phaeo = A.TAXO(:,12);
xx = diatom./counts;
yy = DATA(:,strcmp('dms',headDATA));
iz = strcmp('depth',headDATA);
zr = DATA(:,iz)<20; % depth range
hold on
scatter(xx(zr),yy(zr),50,DATA((zr),iz),'fill'); colormap(cmap);
title(sprintf('r = %0.2f',corr(xx,yy,'rows','pairwise')),'fontsize',16)
xlabel('','fontsize',16)
ylabel('','fontsize',16)
c = colorbar;
ylabel(c,headDATA{iz},'fontsize',16)


%% Correlation between Phaeocystis and diatom counts vs. pigments
ic = 19:2:size(DATA,2);
XC = DATA(:,ic);
removeX = sum(isnan(XC))>30;
XC(:,removeX) = [];
ic(removeX) = [];
iy1 = 11;
iy2 = 2;
R1 = corr(XC,A.TAXO(:,11),'rows','pairwise','type','Spearman');
R2 = corr(XC,A.TAXO(:,2),'rows','pairwise','type','Spearman');
hh = headTAXO(iTAXO);
figure(16),
bar([R1 R2],'barwidth',0.5);
xlim([0 length(ic)+1])
set(gca,'tickdir','out','ticklength',[.005 .005],'xtick',1:length(ic),'xticklabel',[])
for j = 1:length(ic)
    text(j,-0.12,cellstr(headDATA(ic(j))),'rotation',-45,'interpreter','none','color','m')
end
title('Correlation coefficient: Phaeo vs. Pigments','fontsize',16)
legend(hh([iy1 iy2]),'interpreter','none')

%% Quick scatterplot between DMS/Chl and HEX/Chl
cmap = flip(parula);
figure(15),
hex = DATA(:,strcmp(headDATA,'hex'));
tchla = DATA(:,strcmp(headDATA,'tchla'));
xx = hex./tchla;
yy = DATA(:,strcmp('dms',headDATA))./tchla;
iz = strcmp('depth',headDATA);
zr = DATA(:,iz)<30; % depth range
hold on
scatter(xx(zr),yy(zr),50,DATA((zr),iz),'fill'); colormap(cmap);
title(sprintf('r = %0.2f',corr(xx(zr),yy(zr),'rows','pairwise','type','Spearman')),'fontsize',16)
xlabel('','fontsize',16)
ylabel('','fontsize',16)
c = colorbar;
ylabel(c,headDATA{iz},'fontsize',16)


