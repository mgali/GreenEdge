% solar zenith angle and azimuth angle
% https://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF
% hourUTC in decimal units, 1 = 24h

function [sza,saa] = solar_sza_saa(year,doyUTC,hourUTC,timezone,lat,lon)

if mod(year,4), nbdays = 365; else nbdays = 366; end

fractional_year = (2*pi/nbdays).*(doyUTC - 1 + (hourUTC-12)/24);

eqtime = 229.18*(0.000075 + 0.001868 * cos(fractional_year) - 0.032077 * sin(fractional_year)...
    - 0.014615 * cos(2 * fractional_year) - 0.040849 * sin(2 * fractional_year));

decl = 0.006918 - 0.399912*cos(fractional_year) + 0.070257*sin(fractional_year) - 0.006758*cos(2*fractional_year)...
    + 0.000907*sin(2*fractional_year) - 0.002697*cos(3*fractional_year) + 0.00148*sin(3*fractional_year);

time_offset = eqtime + 4*lon - 60*timezone;

tst = hourUTC*24*60 + time_offset;

ha = (tst/4) - 180;

cos_sza = sin(lat).*sin(decl) + cos(lat).*cos(decl).*cos(ha);
sza = acos(cos_sza); % zenith angle

cos_180_minus_aa = (sin(lat).*cos(sza) - sin(decl))./(cos(lat).*sin(sza));
saa = 180 - acos(cos_180_minus_aa);
