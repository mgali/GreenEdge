% READ MODIS ATMOSPHERE DATA AND RUN FORTRAN CODE get_ed0 TO PRODUCE DAILY
% IRRADIANCE FILES AT DESIRED TEMPORAL RESOLUTION
% get_ed READS AND INTERPOLATES 5-DIMENSION LUT TABLES CREATED WITH SBDART BY
% SIMON BELANGER WITH DIMENSIONS: 83 WAVELENGTH, 19 SZA, 8 COT, 10 TO3, 7 SURF ALBEDO
% Marti Gali Tapias 2017-05-08
tic
% ------------------------------------- EDIT -------------------------------------

% Define DOY decrements to output Ed for previous days (comment either case 0 or 1)

% Case 0: normal configuration, only daily data searched
stnlist = dlmread('~/Desktop/GreenEdge/Irradiance/input.noclim.IceCamp20152016.txt');
doydecr = 0;

% % Case 1: GreenEdge cruise, backward search at all stations
% stnlist = dlmread('~/Desktop/GreenEdge/Irradiance/input.noclim.GE2016.txt');
% doydecr = [-2 -1 0];

% Where are files?
% mainpath = '/Volumes/taku-njall/MODISA/L3BIN';
mainpath = '~/Desktop/GreenEdge/Irradiance/Ed_RT_MODISA'; % for tests on MBP

% Above (1) or below (0) water surface?
above = 1;

% What output frequency do we want?
hperiod = 3; % integer

% Do we want to interpolate MODIS data in 2D?
interpol = 0;

% --------------------------------- END EDIT -------------------------------------

if above
    % Above surface: Ed0plus LUT
    lutpath = '/Users/martigalitapias/Documents/SBDART/Ed0plus_LUT_5nm_v2.dat';
    outpath = '~/Desktop/GreenEdge/Irradiance/Ed0plus_MODISA_LUT_SurfAlb';
else
    % Below surface: Ed0moins LUT
    lutpath = '/Users/martigalitapias/Documents/SBDART/Ed0moins_LUT_5nm_v2.dat';
    outpath = '~/Desktop/GreenEdge/Irradiance/Ed0moins_MODISA_LUT_SurfAlb';
end

hours = (0+hperiod/2):hperiod:24;
codepath = '~/Desktop/GreenEdge/Irradiance/fortran_src_MGT';
grid.lat = -(-89.5:1:89.5);
grid.lon = -179.5:1:179.5;
grid.Mlon = meshgrid(grid.lon,grid.lat);
grid.Mlat = (meshgrid(grid.lat,grid.lon))';

% ------------------------- UNCOMMENT FOR TESTS ---------------------------
%
stnlist = [69.32548633 -60.96662 10 6 2016 280]; % for tests on MBP
outpath = sprintf('%s/%04i/%03i/interpol_%01i/',mainpath,stnlist(5),stnlist(6),interpol);
%
% -------------------------------------------------------------------------

% MATCH ICE DATA AND/OR ALBEDO FOR GREENEDGE CRUISE AND ICE CAMPS!

% Loop for stations, days searched at each station, hours within each day
for is = 1:size(stnlist,1)
    
    stn.LAT = stnlist(is,1);
    stn.LON = stnlist(is,2);
    stn.Y = stnlist(is,5);
    stn.DOY = stnlist(is,6);
    
    for dd = 1:length(doydecr)
        
        searchY = stn.Y;
        searchDOY = stn.DOY + doydecr(dd);
        if (searchDOY < 1)
            % Case where changing day changes year
            searchY = stn.Y - 1;
            if ~mod(searchY,4)
                % Leap year
                searchDOY = 366;
            else
                % Normal year
                searchDOY = 365;
            end
        end
        
        % Preallocate output
        OUT = nan(83,length(hours));
        
        % Build outfile path
        if above
            outfilepath = sprintf('%s/Ed0plus_%04i_%03i_%0.3f_%0.3f.txt', outpath, searchY, searchDOY, stn.LAT, stn.LON);
        else
            outfilepath = sprintf('%s/Ed0moins_%04i_%03i_%0.3f_%0.3f.txt', outpath, searchY, searchDOY, stn.LAT, stn.LON);
        end
        
        % Build file path
        genericname = sprintf('MYD08_D3.A%04i%03i.051.*.hdf', searchY, searchDOY);
        searchpath = sprintf('%s/%04i/%03i/%s', mainpath, searchY, searchDOY, genericname);
        
        % Find the complete filename
        getfilepath = sprintf('ls %s',searchpath);
        [status,cmdout] = system(getfilepath);
        
        if status
            
            sprintf('File %s not found',searchpath);
            system(sprintf('echo %s >> list_not_found.txt',searchpath));
            
        else
            
            filepath = strcat(cmdout);
            
            % Read variables
            COT = hdfread(filepath,'/mod08/Data Fields/Cloud_Optical_Thickness_Combined_Mean','Index',{[1 1],[1 1],[180 360]});
            TO3 = hdfread(filepath,'/mod08/Data Fields/Total_Ozone_Mean','Index',{[1 1],[1 1],[180 360]});
            CF = hdfread(filepath,'/mod08/Data Fields/Cloud_Fraction_Day_Mean','Index',{[1 1],[1 1],[180 360]});
            
            % Either interpolate MODIS data in 2D or find direct pixel match
            if interpol
                stn.COT = interp2(grid.Mlon,grid.Mlat,COT,stn.LON,stn.LAT);
                stn.TO3 = interp2(grid.Mlon,grid.Mlat,TO3,stn.LON,stn.LAT);
                stn.CF = interp2(grid.Mlon,grid.Mlat,CF,stn.LON,stn.LAT);
            else
                ilat = 90 - floor(stn.LAT);
                if ilat == 0, ilat = 1; end
                ilon = 180 + ceil(stn.LON);
                if ilon == 0, ilon = 1; end
                stn.COT = COT(ilat,ilon);
                stn.TO3 = TO3(ilat,ilon);
                stn.CF = CF(ilat,ilon);
            end
            
            % Convert to right format and units, replace missing data by NaN
            if stn.COT == -9999, stn.COT = nan; end
            stn.COT = double(stn.COT)*0.01;
            if stn.TO3 == -9999, stn.TO3 = nan; end
            stn.TO3 = double(stn.TO3)*0.1;
            if stn.CF == -9999, stn.CF = nan; end
            stn.CF = double(stn.CF)*0.0001;
            
            % Define albedo
            stn.ALB = 0.05;
            
            % Call get_ed0 in hours loop. For some reason the
            % output is called fort.6
            ! rm fort.6
            for ih = 1:length(hours)
                cmd = sprintf('./get_ed0 %i %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f',...
                    above,searchDOY,hours(ih),stn.LAT,stn.LON,stn.COT,stn.TO3,stn.CF,stn.ALB);
                [status,edout] = system(cmd);
                if status
                    error('get_ed0 failed')
                else
                    READ = dlmread(sprintf('%s/fort.6',codepath));
                    READ = READ(:,(end-82):end);
                    OUT(:,ih) = READ';
                end
                [is searchY searchDOY hours(ih) toc]
            end
            
            %             % Option 2 (almost as fast): write an executable each time.
            %             ! rm fort.6
            %             for ih = 1:length(hours)
            %                 ! rm tmp_get_ed0.sh
            %                 fid = fopen('tmp_get_ed0.sh','w');
            %                 line1 = '#!/bin/bash';
            %                 line2 = sprintf('./get_ed0 %i %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f',...
            %                     above,searchDOY,hours(ih),stn.LAT,stn.LON,stn.COT,stn.TO3,stn.CF,stn.ALB);
            %                 fprintf(fid,'%s\n%s\n',line1,line2);
            %                 ! chmod 755 tmp_get_ed0.sh
            %                 [status,edout] = system('./tmp_get_ed0.sh');
            %                 if status
            %                    error('get_ed0 failed')
            %                 else
            %                     READ = dlmread(sprintf('%s/fort.6',codepath));
            %                     READ = READ(:,(end-82):end);
            %                     OUT(:,ih) = READ';
            %                 end
            %                 [is searchY searchDOY hours(ih) toc]
            %             end
        end
        
        % Write output
        dlmwrite(outfilepath,OUT,'delimiter',' ','precision','%2.4f');
        
    end
end
