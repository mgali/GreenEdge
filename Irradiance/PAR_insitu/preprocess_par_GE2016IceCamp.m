% IMPORT PAR DATA FROM GE ICE CAMP 2016, CALCULATE DAILY PAR AND ALSO PAR
% OF THE 2 DAYS PRIOR TO SAMPLING
clc, clear


A = xlsread('~/Desktop/GreenEdge/Irradiance/PAR_insitu/IC2016/PAR_Data_GE2016-ICECAMP_cyber.xlsx','A3:B1683');


% It seems matlab gets the date with a one-day and 1900 years delay, so
% let's correct it
dummy = datevec(A(:,1));
dummy(:,1) = dummy(:,1) + 1900; % add 1900 years to get 2016
dummy2 = datenum(dummy);
mtimeUTC = dummy2 - 1;

% Calculate local time for Qik (is UTC time minus 4h)
mtimeLOC = mtimeUTC - 4/24;

% Plot to check
figure(),
x = datevec(mtimeUTC);
x = yearday(x(:,3),x(:,2),x(:,1),x(:,4),x(:,5),x(:,6));
plot(x,A(:,2),'-k')
xlabel('Day of year (2016)'), ylabel('PAR, µmol photons m^{-2} s^{-1}')
xlim([118 192])

%% Make histogram
figure(),
hist(A(:,2),0:100:2000)
set(gca,'xtick',0:100:2000,'xticklabel',0:100:2000)
xlim([-100 2100])

% % Save data in mat format
% par_ts.mtimeLOC = mtimeLOC;
% par_ts.mtimeUTC = mtimeUTC;
% par_ts.data = A(:,2);
% parunits = 'micromol m-2 s-1';
% save('PAR_Data_GE2016-ICECAMP_cyber.mat','par_ts','parunits')
