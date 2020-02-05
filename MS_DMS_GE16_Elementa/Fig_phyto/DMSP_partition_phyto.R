# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(RColorBrewer)
library(classInt) # for function classIntervals
# Check http://geog.uoregon.edu/GeogR/topics/multiplots01.html

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
prof <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)

# Exporting image?
exportimg <- T
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_phyto/"


# ---------------------
# Data wrangling

# Remove unnecessary DMS variables
prof$dms <- prof$dms_consens_cf68
toinclude <- names(prof)[grep("dms",names(prof), invert = T)]
toinclude <- c(toinclude,"dms","dmspt")
prof <- prof[,toinclude]

# Remove duplicated variables
prof <- prof[,grep("NA",names(prof), invert = T)]

# Remove data where no phyto counts available
prof <- prof[!is.na(prof$Phaeo) & !is.na(prof$dmspt) & !is.na(prof$depth),]


# If phyto counts duplicated, remove duplicates
# prof <- prof[!duplicated(prof$Phaeo),]

# Add surface and SCM categories
prof$scm <- 'SCM'
prof$scm[prof$depth <= 10] <- 'surface'


# ---------------------
# Convert Phaeocystis counts to DMSP to estimate contribution

# Using DMSP cell quota (q)
qdmsp.phaeo <- 10                                                 # fmol DMSP / cell
prof$dmspt_phaeo.counts <- prof$Phaeo * qdmsp.phaeo * 1e-6        # DMSPt-Phaeo
prof$dmspt_phaeo.counts.pc <- 100*prof$dmspt_phaeo.counts/prof$dmspt
print(summary(prof$dmspt_phaeo.counts.pc))
View(prof[,c("stn","depth","dms","dmspt","dmspt_phaeo.counts","dmspt_phaeo.counts.pc","Phaeo")])

# Using Phaeocystis carbon from iFCB
rdmsp.phaeo <- 0.05
prof$dmspt_phaeo.Cifcb <- prof$phaeo_mg_L * rdmsp.phaeo * (1/12) * (1/5) * 1e6
prof$dmspt_phaeo.Cifcb.pc <- 100*prof$dmspt_phaeo.Cifcb/prof$dmspt
print(summary(prof$dmspt_phaeo.Cifcb.pc))
View(prof[,c("stn","depth","dms","dmspt","dmspt_phaeo.Cifcb","dmspt_phaeo.Cifcb.pc","Phaeo","phaeo_mg_L")])

# ---------------------
# Correlations

# rs.sur <- cor(prof[[pg]][prof$scm=='surface'],prof$dmspt[prof$scm=='surface'], use = "pairwise", method = "spearman")
# rs.scm <- cor(prof[[pg]][prof$scm!='surface'],prof$dmspt[prof$scm!='surface'], use = "pairwise", method = "spearman")
# ps.sur <- cor.test(prof[[pg]][prof$scm=='surface'],prof$dmspt[prof$scm=='surface'], use = "pairwise", method = "spearman")
# ps.scm <- cor.test(prof[[pg]][prof$scm!='surface'],prof$dmspt[prof$scm!='surface'], use = "pairwise", method = "spearman")


# ---------------------
# Plot
# Add lines of C_DMSP:C_tot?

# if (exportimg) {png(filename = paste0(opath,"dms_dmspt_",pg,"_counts.png"), width = 6, height = 7, units = 'cm', pointsize = 8, bg = 'white', res = 600, type = 'cairo')}
# 
# if (exportimg) {dev.off()}
