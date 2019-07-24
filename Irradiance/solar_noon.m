% FIND EXACT TIME OF SOLAR NOON
% https://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF

function snoon = solar_noon(year,doy,hour,lon)
if mod(year,4), nbdays = 365; else nbdays = 366; end
fractional_year = (2*pi/nbdays)*(doy - 1 + (hour-12)/24);
eqtime = 229.18*(0.000075 + 0.001868 * cos(fractional_year) - 0.032077 * sin(fractional_year)...
    - 0.014615 * cos(2 * fractional_year) - 0.040849 * sin(2 * fractional_year));
snoon  = (720 - 4*lon - eqtime)/60;