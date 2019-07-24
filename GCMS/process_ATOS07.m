clear, clc

%% Load and preprocess ATOS07 data

pathA07 = '~/Documents/DADES I CAMPANYES JA PUBLICADES/2010 MarChem - ATOS 2007/CTD-tot';
DATA = csvread(sprintf('%s/CTDok_sofre_complet_4GreenEdge.csv',pathA07),1,0);


%% Process

station = DATA(:,1);
slist = 1:nanmax(station);

% Median cp in 150-200 horizon
cbase = log(100./DATA(DATA(:,4)>=150&DATA(:,4)<=200,9))/0.25; % c 150-200, [m-1]
c.median = nanmedian(cbase);
c.q01 = quantile(cbase,0.01);
c.mean = nanmean(cbase);

DATA_CTD = [];

for k = 1:length(slist) % loop on transect stations
    
    kstn = slist(k);
    tmp = DATA(station==kstn,:);
    
    % Remove if z not surface or DMS or cast number absent
    if sum(isnan(tmp(:,23)))<size(tmp,1)
        
        % Extract CTD variables
        stn = tmp(:,1);
        depth = tmp(:,4);
        temp = tmp(:,5);
        sal = tmp(:,7);
        sigt = tmp(:,15);
        Trans = tmp(:,9);
        O2 = tmp(:,8);
        DMS = tmp(:,23);
        DMS(DMS==0) = nan;
        
        % Calculate cp from Transmittance profile
        cp_out = get_cp(Trans,depth,kstn);
        cp = cp_out.q;
        cpsmooth1 = cp_out.s;
        
        % Append DATA_CTD matrix
        DATA_CTD = [DATA_CTD; [stn depth temp sal sigt O2 cp DMS]];
        
    end
end

%% Save

headDATA_CTD = {'stn' 'depth' 'temp' 'sal' 's-t' 'O2' 'cp' 'DMS'};
save('profile_ATOS07.mat','DATA_CTD','headDATA_CTD');