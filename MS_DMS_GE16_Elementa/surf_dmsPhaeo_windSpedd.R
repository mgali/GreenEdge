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
clco <- ""            # Either NULL (plot all stations), "Arctic" or "Atlantic"
doexploreplot <- F
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_profiles_norm_grouped/"
showtable <- "owd_class"    # either empty, owd_class or sic_class

# ---------------------
pal <- colorRampPalette(brewer.pal(9, "Spectral"))              # Color palette (Consistent with Matlab plots)
col <- pal(n = 21)[c(21,18,5)]
plotres <- 600                                                  # resolution dpi
plet <- sapply(letters, paste0, ")")                            # plot letters

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
# dd <- (duplicated(pplot[,c("dmspt","cast","depth")]) | duplicated(pplot[,c("dms","cast","depth")])) & (!is.na(pplot$tchla) | !is.na(pplot$cpsmooth1)) # Does not work well, too many repeated DMSPt get in
# dd <- (duplicated(pplot[,c("dmspt","cast","depth")]) | duplicated(pplot[,c("dms","cast","depth")])) # Does not work well, too many repeated DMSPt get in
dd <- duplicated(pplot[,c("dmspt","dms","cast","depth")]) | is.na(pplot$dmspt)
pplot <- pplot[!dd,]
# View(pplot[,c("stn","cast","AW_ArW_clustering_coefficient","depth","dms","dmspt","cpsmooth1","tchla")])

# f_mydiff <- function(x) {y <- as.logical(c(1,diff(x))); y[is.na(y)] <- T; return(y)}
# ddms <- f_mydiff(pplot$dms)
# ddepth <- f_mydiff(pplot$depth)
# dcast <- f_mydiff(pplot$cast)
# pplot <- pplot[ddms&!is.na(pplot$dmspt)&!is.na(pplot$tchla),]

# Hide data from stn 400 (either entire or just surface)
# pplot[pplot$stn<=400,] <- NA
pplot[pplot$stn<=400 & pplot$depth < 5,] <- NA

# Change units of N2 from s-2 to h-1
pplot$N2 <- sqrt(pplot$N2) * 3600

# Add ratios
pplot$cp2tchla <- pplot$cpsmooth1/pplot$tchla                           # Cp/tchla ratio

