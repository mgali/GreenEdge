# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(RColorBrewer)
library(dplyr)

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
pdir <- 'plots_pigments_vs_DMS_zColorScale/'
# prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
surf.all <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)

# Exporting image?
exportimg <- T
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/"

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

# Add MIZ classification
icecrit1 <- 0.15
icecrit2 <- 0.70
ICE <- surf.all[,c("SICm2d","SICm1d","SICday")]
icemin <- apply(ICE, 1, min, na.rm = T) # Min-max SIC criterion
icemax <- apply(ICE, 1, max, na.rm = T) # Min-max SIC criterion
icemean <- apply(ICE, 1, min, na.rm = T) # Mean concentration criterion
surf.all$sic_class <- NA
surf.all$sic_class[icemax<icecrit1] <- "OW"
surf.all$sic_class[icemin>icecrit2] <- "ICE"
surf.all$sic_class[icemax>=icecrit1 & icemin<=icecrit2] <- "MIZ"

# Merge with profiles to get clustering coefficient and SIC classification
pplot <- merge(x = prof.all, y = surf.all, all.x = T, all.y = F, by = 'stn', suffixes = "")

# Remove duplicated columns with NA in their names
pplot <- pplot[,grep("NA",names(pplot), invert = T)]

# Remove duplicated rows
dd <- duplicated(pplot[,c("dmspt","dms","cast","depth")]) | is.na(pplot$dmspt)
pplot <- pplot[!dd,]

# Hide data from stn 400 (either entire or just surface)
# pplot[pplot$stn<=400,] <- NA
pplot[pplot$stn<=400 & pplot$depth < 5,] <- NA

# Change units of N2 from s-2 to h-1
pplot$N2 <- sqrt(pplot$N2) * 3600

# Select surface data
pplot <- pplot[pplot$depth == 0.7 & !is.na(pplot$Phaeo), c("stn","depth","dms","Phaeo","tchla","cpsmooth1","wsp24","wsc24","fdmsW97c24")]

# Add variables
pplot$dms2phaeo <- pplot$dms/pplot$Phaeo
pplot$kvent <- pplot$fdmsW97c24/pplot$dms

# Print wsp vs dms corr
print(cor.test(pplot$wsp24, pplot$dms, method = "spearman"))

# -------------------------------
# Remove NA
pplot <- pplot[!is.na(pplot$Phaeo),]

ff <- lsfit(pplot$Phaeo, pplot$dms, wt = NULL, intercept = TRUE, tolerance = 1e-07)
# ff <- lsfit(pplot$tchla, pplot$dms, wt = NULL, intercept = TRUE, tolerance = 1e-07)
# ff <- lsfit(pplot$cpsmooth1, pplot$dms, wt = NULL, intercept = TRUE, tolerance = 1e-07)

print(cor.test(pplot$wsp24, ff$residuals, method = "spearman"))
print(cor.test(pplot$kvent, ff$residuals, method = "spearman"))

