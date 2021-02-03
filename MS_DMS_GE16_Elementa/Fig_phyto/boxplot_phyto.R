# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(dplyr)
library(RColorBrewer)

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
surf.all <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)

# Exporting image?
exportimg <- T
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_phyto/"

# Phyto group plotted
pg <- "Phaeo" # Phaeo, diat_cen, diat_pen
xl <- list(diat_pen = "Pennate diatoms (cells/L)",
           diat_cen = "Centric diatoms (cells/L)",
           Phaeo = 'Phaeocystis (cells/L)')

# Rename DMS variable and remove unnecessary DMS variables
prof.all$dms <- prof.all$dms_consens_cf68
toinclude <- names(prof.all)[grep("dms",names(prof.all), invert = T)]
toinclude <- c(toinclude,"dms","dmspt")
prof.all <- prof.all[,toinclude]

# Add clustering coefficient as per station by merging with profiles
pplot <- merge(x = prof.all, y = surf.all, all.x = T, all.y = F, by = 'stn', suffixes = "")

# Remove DMS data points where DMSPt or Phaeocystis are missing, and stations where microscopy counts not done
pplot <- pplot[!is.na(pplot$dmspt) & !is.na(pplot$dms),]
pplot <- pplot[pplot$stn >= 418,]

# Remove duplicated rows
dd <- duplicated(pplot[,c("dmspt","dms","cast","depth")]) & (!is.na(pplot$tchla) & !is.na(pplot$cpsmooth1))
pplot <- pplot[!dd,]

# Add surface and SCM categories
pplot$scm <- 'SCM'
pplot$scm[pplot$depth < 10] <- 'surface'

xl <- list("diat_cen"="Diat_C",
               "diat_pen"="Diat_P",
               "dino_athec"="Dino_A",
               "dino_thec"="Dino_T",
               "chrys"="Chryso",
               "crypt"="Crypto",
               "Phaeo"="Phaeocystis",
               "flag"="Flag_other")
xpos <- (seq(1, length(xbars)))*6.1 - 1.7


if (exportimg) {png(filename = paste0(opath,"FigS1_phytoCounts.png"), width = 16, height = 8, units = 'cm', pointsize = 8, bg = 'white', res = 600, type = 'cairo')}

# ---------------------
# Multipanel setup
m1_2 <- rbind(matrix(data = 1, nrow = 4, ncol = 9), matrix(data = 2, nrow = 4, ncol = 9))
m3_4 <- rbind(matrix(data = 3, nrow = 4, ncol = 5), matrix(data = 4, nrow = 4, ncol = 5))
m <- cbind(m1_2, m3_4, matrix(data = 5, nrow = 8, ncol = 2))
layout(m)
par(oma = c(1,1,0.5,0.5))

# ---------------------
# a) DMSPt vs. Phyto counts correlations
par(mar = c(4,4,1,1))


# ---------------------
# b) DMS vs. Phyto counts correlations



if (exportimg) {dev.off()}

