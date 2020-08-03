# Compute diapycnal DMS transport at the base of the hBD

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
surf.all <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)

# ---------------------
# Rename DMS variable and remove unnecessary DMS variables
prof.all$dms <- prof.all$dms_consens_cf68
toinclude <- names(prof.all)[grep("dms",names(prof.all), invert = T)]
toinclude <- c(toinclude,"dms","dmspt")
prof.all <- prof.all[,toinclude]

# Remove data where no DMS or DMSPt are available
prof.all <- prof.all[(!is.na(prof.all$dms) | !is.na(prof.all$dmspt)) & !is.na(prof.all$depth),]

# !!!!!! Correct DMS and DMSPt in stn 519 surface !!!!!
# prof.all[prof.all$stn==519 & prof.all$depth==0.7,c("dms","dmspt")] <- c(3.93,79.9)
prof.all[prof.all$stn==519 & prof.all$depth==0.7,c("dms","dmspt")] <- c(11.42,79.9)
# ---------------------

# Compute vertical DMS gradient
prof.all$dDMSdz <- c(
  diff(prof.all$dms)/diff(prof.all$depth),
  NA
)
prof.all$dDMSdz[abs(prof.all$dDMSdz)==Inf] <- NA
prof.all$dDMSdz[prof.all$depth>20] <- NA
print(summary(prof.all$dDMSdz))
# View(prof.all[,c("stn","depth","dms","dDMSdz")])

# Compute vertical flux
kzmean <- 5e-5
prof.all$kzDMS <- prof.all$dDMSdz * kzmean *  86400 # µmol m-2 d-1
print(summary(prof.all$kzDMS))

# Compute volumetric rate: need to fill each station's data with hBD: not always meaningful
# Compute specific rate

# ---------------------
# Find highest dDMSdz
ff <- which(prof.all$dDMSdz > quantile(prof.all$dDMSdz, 0.95, na.rm = T))
print(prof.all$dDMSdz[ff])

# Max gradient at stn. 512 between 0.7 and 15 m, mean 4.3 µmol m-4
# hBD is 11.3 m, MLD0.1 is 5 m, MLD0.03 is 1.6 m
# DMS is 12 at 0.7 m, 36 at 5 m, 74 at 15 m
kzDMS <- 4.3 * 2e-5 * 86400 # µmol m-2 d-1
print(paste0("Vertical transport = ",round(kzDMS,3)," µmol m-2 d-1")) # Should not be higher than the sea-air flux! Estimated at 12 µmol m-2 d-1

# Volumetric rate
hbd <- 11.3
vkzDMS <- kzDMS / hbd # µmol m-3 d-1
print(paste0("Volumetric vertical transport = ",round(vkzDMS,3)," µmol m-3 d-1")) # Should not be higher than the sea-air flux! Estimated at 12 µmol m-2 d-1

# Specific rate
meanDMS <- 0.7 + 4.3*5
print(paste0("Rate constant = ",round(vkzDMS/meanDMS,3)))

