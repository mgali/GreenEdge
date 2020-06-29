# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(dplyr)
library(MASS)
library(vegan)
install.packages("vegan")

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
prof <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)

# Output settings
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/RDA/"
exportimg <- T

# Rename DMS variable and remove unnecessary DMS variables
prof$dms <- prof$dms_consens_cf68
toinclude <- names(prof)[grep("dms",names(prof), invert = T)]
toinclude <- c(toinclude,"dms","dmspt")
prof <- prof[,toinclude]
prof <- prof[,names(prof)[grep("cf",names(prof), invert = T)]]

# Remove data where no DMS or DMSPt are available (needed for assignment of new values in next step)
prof <- prof[(!is.na(prof$dms) | !is.na(prof$dmspt)) & !is.na(prof$depth),]

# !!!!!! Correct DMS and DMSPt in stn 519 surface !!!!!
prof[prof$stn==519 & prof$depth==0.7,c("dms","dmspt")] <- c(11.42,79.9)

# # Count proportion of empty cells in each column (mostly used to exclude some pigments or group them into functional units as DD or VAZ xanthophylls)
# noNAcount <- sapply(prof, function(x) round( sum(!is.na(x) & (!is.na(prof$dms) | !is.na(prof$dmspt))), digits = 2) )

# Group DD and VAZ cycles
prof$dd <- rowSums(prof[,c("diadino","diato")], na.rm = T)
prof$vaz <- rowSums(prof[,c("zea","anthera","viola")], na.rm = T)

# ------------------------------------------------------------------------
# Prepare regression and run in loop
# Pigments with less than 50% available data

yvarS <- c("dmspt","dms")
xvarS <- list(physics = c("temp","sal","N2","cpsmooth1","anp","par_d_p24h_ein_m_2_day_1","dmspt"),
              pigments = c("chlc3","chlc2group","chldaSUM","peri","phdaSUM","but","fuco","neo","pras.1","hex",
                           "dd","allo","lut","chlb","tchla","phytnaSUM","tcar","but19_like","vaz"))
