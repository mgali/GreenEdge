% ASSIGN ICE COVER AND ALBEDO CATEGORY
%
% Marti Gali Tapias 2017-05-26, v2


function [type,ialbedo,talbedo] = f_get_albedo_v3(doy,sic,snow)

% Parameter definitions
calbedo = 0.85; % cold snow albedo
malbedo = 0.7; % slope of albedo vs. sic for melting ice/snow
walbedo = 0.1; % open water albedo, intercept of albedo vs. sic linear relationship
crit_sic1 = 0.85; % sea ice concentration threshold, 1
crit_snow = 0.1; % snow depth threshold in m
crit_doy1 = 152; % date when melting-ponding begins

% Preallocate
type = nan(size(doy)); % ice/snow cover type
ialbedo = nan(size(doy)); % ice/snow albedo
talbedo = nan(size(doy)); % total albedo

% If sic >= crit_sic1 and either snow is present (>=0.1 m) % or doy is 
% earlier than crit_doy1 (June 1st), assign type1, albedo is that 
% of cold snow according to Perovich 2012.
coldsnow = doy < crit_doy1 & (sic >= crit_sic1 | (~isnan(snow) & snow >= crit_snow));
type(coldsnow) = 1;
ialbedo(coldsnow) = calbedo;
talbedo(coldsnow) = ialbedo(coldsnow).*sic(coldsnow) + walbedo*(1 - sic(coldsnow));

% If not coldsnow: albedo is linear function of sic
melting = ~coldsnow;
type(melting) = 2;
ialbedo(melting) = malbedo;
talbedo(melting) = ialbedo(melting).*sic(melting) + walbedo*(1 - sic(melting));

% If not assigned
unassigned = ~coldsnow & ~melting;
type(unassigned) = 0;