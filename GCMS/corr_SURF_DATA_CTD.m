clear, clc

load('surf_profile_ALL_4plot.mat','SURF','headSURF','DATA_CTD','headDATA_CTD')


%% --------------- Explore relationships in DATA_CTD matrix ---------------

%% Linear and log space correlations

% % For DATA_CTD matrix (profiles)
% include = {'dms_consens_cf68' 'depth' 'temp' 'sal' 'sigt' 'N2' 'O2' 'CDOM' 'NO3' 'cpsmooth1'};
% shortname = {'DMS' 'Z' 'T' 'S' 's-t' 'N2' 'O2' 'CDOM' 'NO3' 'cp'};
% Mall = [];
% XCORR = DATA_CTD;
% headX = headDATA_CTD;

% SURF matrix
include = headSURF([7:12 14:16 18:23 25 26]);
shortname = {'dms' 'SIC2' 'SIC1' 'SICd' 'WS72' 'WS24' 'SST' 'SAL' 'Fdms' 's-t' 'O2' 'CDOM' 'NO3' 'cp' 'XLD' 'N2m' 'zN2m'};
Mall = [];
XCORR = SURF;
headX = headSURF;

for j = 1:length(include)
    ivar = include{j};
    Mall = [Mall XCORR(:,strcmp(headX,ivar))];
end
[R.a,P.a] = corr(Mall,'rows','pairwise');
[R.at,P.at] = corr(real(log10(Mall)),'rows','pairwise');

%% Check histograms
% for j = 1:length(include)
%     figure(j),
%     hist(Mall(:,j),10)
%     xlabel(include{j},'interpreter','none')
%     print(j,sprintf('hist_%s.png',include{j}),'-dpng','-r300')
%     close(j)
% end

%% Plot correlation matrices: lin and log side by side
cmap = flip(brewermap(21,'RdBu'));
cmap(9:11,:) = ones(3,3);

X = tril(R.a);
Y = tril(P.a);
Y(X==1) = 1;
X(X==1) = 0;

figure(1), clf
set(gcf,'units','centimeters','position',[2 2 36 16])

subplot(121)
imagesc(X)
colormap(cmap)
set(gca,'xtick',1:length(include),'xticklabel',shortname,'ytick',1:length(include),'yticklabel',shortname)
for j = 1:length(include)
    for k = 1:length(include)
        if X(j,k) ~= 0 && X(j,k) ~= 1
            if Y(j,k) <= 0.05 && Y(j,k) > 0.01
                text(k-0.3,j,sprintf('%0.2f',X(j,k)),'color','w','fontsize',10)
            elseif Y(j,k) <= 0.01
                text(k-0.32,j,sprintf('%0.2f',X(j,k)),'fontweight','bold','color','w','fontsize',11)
            end
        end
    end
end
title(sprintf('Linear space correlations \n p < 0.05 only'),'fontsize',16)


X = tril(R.at);
Y = tril(P.at);
Y(X==1) = 1;
X(X==1) = 0;

subplot(122)
imagesc(X)
colormap(cmap)
set(gca,'xtick',1:length(include),'xticklabel',shortname,'ytick',1:length(include),'yticklabel',shortname)
for j = 1:length(include)
    for k = 1:length(include)
        if X(j,k) ~= 0 && X(j,k) ~= 1
            if Y(j,k) <= 0.05 && Y(j,k) > 0.01
                text(k-0.3,j,sprintf('%0.2f',X(j,k)),'color','w','fontsize',10)
            elseif Y(j,k) <= 0.01
                text(k-0.32,j,sprintf('%0.2f',X(j,k)),'fontweight','bold','color','w','fontsize',11)
            end
        end
    end
end
title(sprintf('Log space correlations \n p < 0.05 only'),'fontsize',16)

%% Plot correlation matrix
cmap = flip(brewermap(21,'RdBu'));
cmap(9:11,:) = ones(3,3);

X = tril(R.a);
Y = tril(P.a);
Y(X==1) = 1;
X(X==1) = 0;

figure(11), clf
set(gcf,'units','centimeters','position',[2 2 20 20])

imagesc(X)
colormap(cmap)
set(gca,'xtick',1:length(include),'xticklabel',shortname,'ytick',1:length(include),'yticklabel',shortname)
for j = 1:length(include)
    for k = 1:length(include)
        if X(j,k) ~= 0 && X(j,k) ~= 1
            if Y(j,k) <= 0.05 && Y(j,k) > 0.01
                text(k-0.3,j,sprintf('%0.2f',X(j,k)),'color','w','fontsize',10)
            elseif Y(j,k) <= 0.01
                text(k-0.32,j,sprintf('%0.2f',X(j,k)),'fontweight','bold','color','w','fontsize',11)
            end
        end
    end
end
title(sprintf('Linear space correlations \n p < 0.05 only'),'fontsize',16)
