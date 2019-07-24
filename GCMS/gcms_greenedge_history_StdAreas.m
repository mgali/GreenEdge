% Analyze correlation between DMS 62 and peak areas of standards 65-68
% Eventually correct 62 areas back in xls spreadsheet
clear, clc, close all

%% Load and pre-process
load history_areas.mat

Aorig = xlsread('~/Desktop/GreenEdge/GCMS/gcms_greenedge_Marti.xlsx','Niskin_underway','A2:AD135');
Aorig(:,16:18) = Aorig(:,16:18)*20./(Aorig(:,14)*ones(1,3)); % Correct for injection volume

%% 1. Plot areas of 65 and 68 standards over time
date = datenum(Aorig(:,1:3));
% /////////////////
% date = date + Aorig(:,29); % uncomment if wanting to sort by analysis date
% /////////////////
DUMMY = [Aorig date];
DUMMY = sortrows(DUMMY,size(DUMMY,2));
Asort = DUMMY(:,1:end-1);
date = DUMMY(:,end);

%%
figure(), hold on
col = 18; % choose 16, 17 or 18, black, blue or green points, m/z 62, 65 or 68
if col == 16, ticky = 0:2e5:14e5; dotcolor = 'k'; textposy = -200000; laby = 62;
elseif col == 17, ticky = 0:5e3:3.5e4; dotcolor = 'b'; textposy = -5000; laby = 65;
elseif col == 18, ticky = 0:5e4:3e5; dotcolor = 'g'; textposy = -40000; laby = 68;
end
yplot = Asort(:,col);
dateplot = date;
dateplot(isnan(yplot)) = [];
diffdate = diff(dateplot);
dateax = [dateplot(diffdate~=0); dateplot(end)];
dateaxstr = cellstr(datestr(dateax));
yplot(isnan(yplot)) = [];

% Axis for individual data points with jitter
xplotjitter = 1;
c = 1;
for j = 1:length(dateplot)-1
    if ~diffdate(j)
    xplotjitter = [xplotjitter; xplotjitter(end)];
    elseif diffdate(j)
        xplotjitter = [xplotjitter; xplotjitter(end)+1];
    end
end
xplotjitter = xplotjitter + randn(length(xplotjitter),1)/10;

% Means, medians
plot(1:16,nanmean(yplot)*ones(1,16),'-r','linewidth',1)
plot(1:16,nanmedian(yplot)*ones(1,16),'--r','linewidth',1)

plot(xplotjitter,yplot,'.','markersize',15,'color',dotcolor), hold on
legend('Mean','Median','All data')
boxplot(yplot,dateplot,'colors',[.5 .5 .5],'plotstyle','traditional',...
    'outliersize',0.001,'whisker',1.5,'widths',[.6 .6 .6 .6],'notch','off');
set(gca,'xtick',1:length(dateaxstr),'xticklabel',cell(15,1),'ytick',ticky,'ylim',[min(ticky) max(ticky)])
for j = 1:length(dateaxstr)
   text(j-1.3,textposy,dateaxstr{j},'rotation',45)
end

ylabel(sprintf('m/z %0.0f area',laby),'fontsize',14)



%% 2. Regress 68 and 65 vs 62 for the relevant range (that with maximal correlation)
% as found in section 3.

% ///// Comment if not correcting for 62-driven variability ///// 
% Correct m/z 68
cut = 400000;
useinreg = Asort(:,16) <= cut;
lmean = nanmean(Asort(useinreg,18));
y = Asort(useinreg,18)-lmean;
X = [ones(sum(useinreg),1) Asort(useinreg,16)];
[b,bi,r,ri,st] = regress(y,X); b, (bi(:,2)-bi(:,1))/2, st([1 3])
correction = X*b;
Asort(useinreg,18) = Asort(useinreg,18) - correction;

% Correct m/z 65
cut = 500000;
useinreg = Asort(:,16) <= cut;
lmean = nanmean(Asort(useinreg,17));
y = Asort(useinreg,17)-lmean;
X = [ones(sum(useinreg),1) Asort(useinreg,16)];
[b,bi,r,ri,st] = regress(y,X); b, (bi(:,2)-bi(:,1))/2, st([1 3])
correction = X*b;
Asort(useinreg,17) = Asort(useinreg,17) - correction;
% ////////////////////////////////////////////////////////////////////// 


%% 3. Scatterplot between different peak areas. Correlation between m/z 68 or
% 65 and 62 peak areas with increasing maximum cutoff 62 PA.

max62PA = (2.5e4:2.5e4:14e5)';
corr6562 = nan(length(max62PA),2);
corr6862 = nan(length(max62PA),2);

max68PA = (4e4:2e3:3e5)';
corr6568 = nan(length(max68PA),2);

xcorr = nansum(Asort(:,16),2); % used to correlate sum of 3 masses to each iso

for ii = 1:length(max62PA)
    select = xcorr <= max62PA(ii);
    if sum(select)>0,
        [r,p] = corr(xcorr(select==1),Asort(select==1,5),'rows','pairwise'); r(p>0.05)=nan;
        corr6562(ii,:) = [r,p];
        [r,p] = corr(xcorr(select==1),Asort(select==1,6),'rows','pairwise'); r(p>0.05)=nan;
        corr6862(ii,:) = [r,p];
    end
end
for ii = 1:length(max68PA)
    select = Asort(:,6) <= max68PA(ii);
    if sum(select)>0,
        [r,p] = corr(Asort(select==1,5),Asort(select==1,6),'rows','pairwise'); r(p>0.05)=nan;
        corr6568(ii,:) = [r,p];
    end
end

figure(65), hold on
scatter(xcorr,Asort(:,17),'filled','c');
ax1 = gca;
plot(max62PA,corr6562*3.5e4,'-','linewidth',3);
ax2 = axes('Position',get(ax1,'Position'),...
    'XTickLabel',[],... %'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k');
linkaxes([ax1 ax2],'x');
xlabel('m/z 62 area','fontsize',14)
ylabel(ax1,'m/z 65 area','fontsize',14)
ylabel(ax2,'Correlation, p','fontsize',14)


figure(68), hold on
scatter(xcorr,Asort(:,18),'filled','g')
ax1 = gca;
plot(max62PA,corr6862*3e5,'-','linewidth',3)
ax2 = axes('Position',get(ax1,'Position'),...
    'XTickLabel',[],... %'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k');
linkaxes([ax1 ax2],'x');
xlabel('m/z 62 area','fontsize',14),
ylabel(ax1,'m/z 68 area','fontsize',14)
ylabel(ax2,'Correlation, p','fontsize',14)


figure(58), hold on
scatter(Asort(:,18),Asort(:,17),'filled','K')
ax1 = gca;
plot(max68PA,corr6568*3.5e4,'-','linewidth',3)
ax2 = axes('Position',get(ax1,'Position'),...
    'XTickLabel',[],... %'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k');
linkaxes([ax1 ax2],'x');
xlabel('m/z 68 area','fontsize',14),
ylabel(ax1,'m/z 65 area','fontsize',14)
ylabel(ax2,'Correlation, p','fontsize',14)


% max62PA(corr6562==nanmax(corr6562(:,1)))
% max62PA(corr6862==nanmax(corr6862(:,1)))

sprintf('Ratio of means = %0.2f, Mean ratio = %0.2f',...
nanmean(Asort(:,18))/nanmean(Asort(:,17)),...
nanmean(Asort(:,18)./Asort(:,17)))


% Mean standard concentrations
Aorig(:,16:18) = Aorig(:,16:18).*(Aorig(:,14)*ones(1,3))/20; % De-correct for injection volume
sprintf('Mean 68 concentration = %0.2f, theoretical 5.0',nanmean(1e9*(20./Asort(:,14)).*10.^(0.968*log10(Asort(:,18))-13.26)))
sprintf('Mean 65 concentration = %0.2f, theoretical 0.54',nanmean(1e9*(20./Asort(:,14)).*10.^(0.968*log10(Asort(:,17))-13.26)))


