# Calculate relative weight of MIZ EDMS compared to open ocean EDMS, in 5 degrees latitude bands

datapath <- '~/Desktop/Artic_DOSES/Arctic_DOSES_timeseries/time_series_3Dmatrix/'

# Read data
ooedms <- list(year = read.csv(paste0(datapath,'INTFLUX_MIZcompare_5degBins_year_TS_fdmsNO_gsm_filled_AM20032016_FDMS_W97_8D_28km.mat.csv'), header = F),
               mjja = read.csv(paste0(datapath,'INTFLUX_MIZcompare_5degBins_MJJA_TS_fdmsNO_gsm_filled_AM20032016_FDMS_W97_8D_28km.mat.csv'), header = F),
               jj = read.csv(paste0(datapath,'INTFLUX_MIZcompare_5degBins_JJ_TS_fdmsNO_gsm_filled_AM20032016_FDMS_W97_8D_28km.mat.csv'), header = F))

mizedms <- list(A = read.csv(paste0(datapath,'MIZ_INTMIZ_caseA_5_x_years_65N85N_20032014.csv'), header = F),
                B = read.csv(paste0(datapath,'MIZ_INTMIZ_caseB_5_x_years_65N85N_20032014_8d.csv'), header = F))

# Scale mizedms by fixedflux, ndays and units conversion factor (INTMIZ contains only the sum of MIZ area in a given period in km2)
fixedflux <- 1 # µmol m-2 s-1 = mol km-2 d-1
ndays <- 8 # days
cf <- 32e-9 # Gg S / mol S

# lapply(mizedms, function(x) x * fixedflux * ndays * cf)

# # Restrict ooedms to same latbands (65-70, 70-75,75-80,80,85) and years
# lapply(ooedms, function(x) {x <- x[5:8,1:12]})

# # Name columns (latbands) and rows (years)
# lapply(ooedms, function(x) rownames(x) <- seq(45,85,5))
# lapply(ooedms, function(x) colnames(x) <- seq(2003,2016,1))


# Quick calculation of mean and std emission from comparable domains

selmiz <- 'B'
selperiod <- 'mjja'

for (j in 1:4) {
  
  mmiz <- mean(as.numeric(mizedms[[selmiz]][j,]), na.rm = T) * fixedflux * ndays * cf
  moo <- mean(as.numeric(ooedms[[selperiod]][j+4,]), na.rm = T)
  
  print(paste0('Mean MIZ EDMS = ',mmiz))
  print(paste0('Mean OO EDMS = ',moo))
  print(paste0('Mean MIZ/OO % = ',100*mmiz/moo))
  
}