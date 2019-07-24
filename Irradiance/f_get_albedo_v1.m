% ASSIGN ICE COVER AND ALBEDO CATEGORY
% This v1 uses extremely simple parameterization with 3 ice categories and
% their corresponding albedos:
%
%   cold snow, albedo = 0.85
%   melting snow = 0.7
%   ponded and melting ice, albedo = 0.5
%   open water, albedo = 0.07
%
% After Perovich et al. 2002 doi:10.1029/2000JC000438, Perovich &
% Polashenski 2012 doi:10.1029/2012GL051432 and Flocco et al. 2012 doi:10.1029/2012JC008195
%
% Marti Gali Tapias 2017-05-17, v1


function [type,ialbedo,talbedo] = f_get_albedo_v1(doy,sic,snow)

% Parameter definitions
calbedo = 0.85; % cold snow albedo
msalbedo = 0.7; % melting snow albedo
mpialbedo = 0.5; % melting/ponded ice albedo (increased 0.45->0.5 after low bias obesrved in v1.1)
walbedo = 0.07; % open water albedo
crit_sic1 = 0.9; % sea ice concentration threshold, 1
crit_sic2 = 0.75; % sea ice concentration threshold, 2
crit_snow = 0.1; % snow depth threshold in m
crit_doy1 = 152; % date when melting-ponding begins
crit_doy2 = 182; % date when ponding starts clearly decelerating at pan-Arctic level

% Preallocate
type = nan(size(doy)); % ice/snow cover type
ialbedo = nan(size(doy)); % ice/snow albedo
talbedo = nan(size(doy)); % total albedo

% If sic >= crit_sic1 and either snow is present (>=0.1 m) % or doy is 
% earlier than crit_doy1 (June 1st), assign type1, albedo is that 
% of cold snow according to Perovich 2012.
coldsnow = sic >= crit_sic1 & (doy < crit_doy1 | (~isnan(snow) & snow >= crit_snow));
type(coldsnow) = 1;
ialbedo(coldsnow) = calbedo;
talbedo(coldsnow) = ialbedo(coldsnow).*sic(coldsnow) + walbedo*(1 - sic(coldsnow));

% If not coldsnow but sic >= crit_sic2 and either we are in June or snow > crit_snow
msnow = ~coldsnow & sic >= crit_sic2 & (doy < crit_doy2 | (~isnan(snow) & snow >= crit_snow));
type(msnow) = 2;
ialbedo(msnow) = msalbedo;
talbedo(msnow) = ialbedo(msnow).*sic(msnow) + walbedo*(1 - sic(msnow));

% If neither coldsnow nor melting snow
mpice = ~isnan(sic) & ~isnan(doy) & ~coldsnow & ~msnow;
type(mpice) = 3;
ialbedo(mpice) = mpialbedo;
talbedo(mpice) = ialbedo(mpice).*sic(mpice) + walbedo*(1 - sic(mpice));

% If not aissgned
unassigned = ~coldsnow & ~msnow & ~mpice;
type(unassigned) = 0;