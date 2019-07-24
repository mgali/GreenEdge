% FIND OPTIMAL CALIBRATION
clear

%% Read data and append

data1 = xlsread('~/Desktop/GreenEdge/GCMS/gcms_greenedge_Marti.xlsx','Niskin_underway','A2:Z135');
data1(:,24) = nan;

data2 = xlsread('~/Desktop/GreenEdge/GCMS/gcms_greenedge_Marti.xlsx','data4MatlabCalibration','A2:Z30');
data2(1:4,:) = nan;

% data = [data1; data2];
data = data2;

%% Assign variable names
injvol = data(:,14); % injection volume in mL
pa62 = data(:,16); % peak areas
pa65 = data(:,17);
pa68 = data(:,18);
ec62 = data(:,24); % expected concentrations in nmol/L
ec65 = data(:,25);
ec68 = data(:,26);

%% Define x and y for fit
x = [pa62; pa65; pa68];
iso = [62*ones(size(pa62)); 65*ones(size(pa65)); 68*ones(size(pa68))];
injvol = 1e-3*repmat(injvol,[3 1]); % in L
ec = 1e-9*[ec62; ec65; ec68]; % in mol/L
y = 1e12*ec.*injvol; % in pmol = 1e12*mol

%% Fit: pmol S vs PA. Calculate stats

xy = [x y];
xy = real(log10(xy));

% Removing outliers has no influence on linear fit
% outliers = [1 2 136 137 138 141 142 143 148 150 151 163 308]; xy(outliers,:) = nan;

[xy,index] = sortrows(xy,1);
iso = iso(index);
iso(isnan(xy(:,1))|isnan(xy(:,2)),:) = [];
xy(isnan(xy(:,1))|isnan(xy(:,2)),:) = [];
x = xy(:,1);
y = xy(:,2);

f = fit(x,y,'poly1'); f
r2 = (corr(f(x),y))^2
p11 = predint(f,x,0.95,'observation','off');
p12 = predint(f,x,0.95,'observation','on');
p21 = predint(f,x,0.95,'functional','off');
p22 = predint(f,x,0.95,'functional','on');

% outliers = find(y<p11(:,1) | y>p11(:,2));

% % Regular linear fit
% [b,bi,r,ri,st] = regress(y,[ones(size(x)) x]); b, bi, st(1)

%% Plot
figure(1), clf,
set(gcf,'units','centimeters','position',[5 5 30 30])

subplot(2,2,1)
plot(f,x,y,'ok'), hold on
legend('location','northwest');
plot(x(iso==62),y(iso==62),'.r','markersize',30)
plot(x(iso==65),y(iso==65),'.b','markersize',10)
plot(x(iso==68),y(iso==68),'.g','markersize',10)
plot(x,p11,'m--')
title('Nonsimultaneous observation bounds','Color','m')

subplot(2,2,2)
plot(f,x,y,'ok'), hold on
legend('location','northwest');
plot(x(iso==62),y(iso==62),'.r','markersize',30)
plot(x(iso==65),y(iso==65),'.b','markersize',10)
plot(x(iso==68),y(iso==68),'.g','markersize',10)
plot(x,p12,'m--'),
title('Simultaneous observation bounds','Color','m')

subplot(2,2,3)
plot(f,x,y,'ok'), hold on
legend('location','northwest');
plot(x(iso==62),y(iso==62),'.r','markersize',30)
plot(x(iso==65),y(iso==65),'.b','markersize',10)
plot(x(iso==68),y(iso==68),'.g','markersize',10)
plot(x,p21,'m--'),
title('Nonsimultaneous functional bounds','Color','m')

subplot(2,2,4)
plot(f,x,y,'ok'), hold on
legend('location','northwest');
plot(x(iso==62),y(iso==62),'.r','markersize',30)
plot(x(iso==65),y(iso==65),'.b','markersize',10)
plot(x(iso==68),y(iso==68),'.g','markersize',10)
plot(x,p22,'m--'),
title('Simultaneous functional bounds','Color','m')
