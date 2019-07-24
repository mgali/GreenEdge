% ASSIGN ICE COVER CATEGORY
% This v0 uses extremely simple parameterization with 4 ice categories and
% their corresponding albedos:
%
%   cold snow, albedo = 0.85
%   melting snow, ponded and melting ice, albedo = 0.5
%   open water, albedo = 0.05
%
% After Perovich et al. 2002 doi:10.1029/2000JC000438, Perovich &
% Polashenski 2012 doi:10.1029/2012GL051432 and Flocco et al. 012 doi:10.1029/2012JC008195
%
% Marti Gali Tapias 2017-05-17, v0


function [type,ialbedo,talbedo] = f_get_albedo_v0(doy,sic,snow)

% Parameter definitions
calbedo = 0.85; % cold snow albedo
malbedo = 0.5; % melting/ponded ice albedo
walbedo = 0.05; % open water albedo
crit_sic = 0.9; % sea ice concentration threshold
crit_snow = 0.1; % snow depth threshold in m
crit_doy = 152; % date when melting-ponding begins

% Preallocate
type = nan(size(doy)); % ice/snow cover type
ialbedo = nan(size(doy)); % ice/snow albedo
talbedo = nan(size(doy)); % total albedo

% If sic >= crit_sic, or if sic < crit_sic and snow is present (>=0.1 m) 
% or doy is earlier than crit_doy (June 1st), assign type1, albedo is that 
% of cold snow according to Perovich 2012.
coldsnow = sic >= crit_sic | doy < crit_doy | (~isnan(snow) & snow >= crit_snow);
type(coldsnow) = 1;
ialbedo(coldsnow) = calbedo;
talbedo(coldsnow) = ialbedo(coldsnow).*sic(coldsnow) + walbedo*(1 - sic(coldsnow));

% If not coldsnow
melting = ~isnan(sic) & ~isnan(doy) & ~coldsnow;
type(melting) = 2;
ialbedo(melting) = malbedo;
talbedo(melting) = ialbedo(melting).*sic(melting) + walbedo*(1 - sic(melting));

% If not aissgned
unassigned = ~melting & ~coldsnow;
type(unassigned) = 0;