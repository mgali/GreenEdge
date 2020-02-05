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

# Function to compute carbon content of phytoplankton
MendenDeuer_Lessard <- function(esd_microm, cell_conc, is_diatom) {
  vol <- (4/3) * pi * (esd_microm/2) # µm^3
  if (is_diatom) {
    pgC <- 0.288 * vol^0.811
  } else {
    pgC <- 0.216 * vol^0.939
  }
  return(pgC * cell_conc)
}

# ---------------------
# Data wrangling

# Remove unnecessary DMS variables
prof$dms <- prof$dms_consens_cf68
toinclude <- names(prof)[grep("dms",names(prof), invert = T)]
toinclude <- c(toinclude,"dms","dmspt")
prof <- prof[,toinclude]

# Remove duplicated variables appended with NA
prof <- prof[,grep("NA",names(prof), invert = T)]

# Remove data where no phyto counts available
prof <- prof[!is.na(prof$Phaeo) & !is.na(prof$dmspt) & !is.na(prof$depth),]
prof <- prof[prof$depth!=0.1 & !is.na(prof$prescreen) ,]

# Add surface and SCM categories
prof$scm <- 'SCM'
prof$scm[prof$depth <= 10] <- 'surface'


# ---------------------
# Convert Phaeocystis and diatom counts to DMSP to estimate contribution

# Using DMSP cell quota (q)
qdmsp.phaeo <- 10                                                 # fmol DMSP / cell
prof$dmspt_phaeo.counts <- prof$Phaeo * qdmsp.phaeo * 1e-6        # DMSPt-Phaeo
prof$dmspt_phaeo.counts.pc <- 100*prof$dmspt_phaeo.counts/prof$dmspt
print(summary(prof$dmspt_phaeo.counts.pc))

# Using Phaeocystis and diatoms carbon
prof$Phaeo_mgC_L.microscopy <- MendenDeuer_Lessard(esd_microm = 5, cell_conc = prof$Phaeo, is_diatom = F) * 1e-9
rdmsp.phaeo <- 0.10
prof$dmspt_phaeo.Cmic <- prof$Phaeo_mgC_L.microscopy * rdmsp.phaeo * (1/12) * (1/5) * 1e6
prof$dmspt_phaeo.Cmic.pc <- 100 * prof$dmspt_phaeo.Cmic / prof$dmspt
print(summary(prof$dmspt_phaeo.Cmic.pc))

prof$diat_cen_mgC_L.microscopy <- MendenDeuer_Lessard(esd_microm = 100, cell_conc = prof$diat_cen, is_diatom = T) * 1e-9
prof$diat_pen_mgC_L.microscopy <- MendenDeuer_Lessard(esd_microm = 20, cell_conc = prof$diat_pen, is_diatom = T) * 1e-9
prof$diat_mgC_L.microscopy <- prof$diat_cen_mgC_L.microscopy + prof$diat_pen_mgC_L.microscopy
rdmsp.diat <- 0.02
prof$dmspt_diat.Cmic <- prof$diat_mgC_L.microscopy * rdmsp.diat * (1/12) * (1/5) * 1e6
prof$dmspt_diat.Cmic.pc <- 100 * prof$dmspt_diat.Cmic / prof$dmspt
print(summary(prof$dmspt_diat.Cmic.pc))

# # Using Phaeocystis carbon from iFCB: WEIRD
# rdmsp.phaeo <- 0.05
# prof$dmspt_phaeo.Cifcb <- prof$phaeo_mg_L * rdmsp.phaeo * (1/12) * (1/5) * 1e6
# prof$dmspt_phaeo.Cifcb.pc <- 100*prof$dmspt_phaeo.Cifcb/prof$dmspt
# print(summary(prof$dmspt_phaeo.Cifcb.pc))


# OUTPUT
OUT <- prof[,c("stn","depth",
               "Phaeo","diat_cen","diat_pen",
               "Phaeo_mgC_L.microscopy","diat_cen_mgC_L.microscopy","diat_pen_mgC_L.microscopy","diat_mgC_L.microscopy",
               "dmspt_phaeo.counts","dmspt_phaeo.Cmic","dmspt_diat.Cmic",
               "dmspt",
               "dmspt_phaeo.counts.pc","dmspt_phaeo.Cmic.pc","dmspt_diat.Cmic.pc")]
View(OUT)
write.csv(x = OUT, file = paste0(opath,"Fraction_DMSPt_Phaeo_diat.csv"), row.names = F)

# ---------------------
# Correlations

Corr <- list(
  sur = list(
    rs = cor(prof[[pg]][prof$scm=='surface'],prof$dmspt[prof$scm=='surface'], use = "pairwise", method = "spearman"),
    ps = cor.test(prof[[pg]][prof$scm=='surface'],prof$dmspt[prof$scm=='surface'], use = "pairwise", method = "spearman")
  ),
  scm = list(
    rs = cor(prof[[pg]][prof$scm!='surface'],prof$dmspt[prof$scm!='surface'], use = "pairwise", method = "spearman"),
    ps = cor.test(prof[[pg]][prof$scm!='surface'],prof$dmspt[prof$scm!='surface'], use = "pairwise", method = "spearman")
  )
)


# ---------------------
# Plot
# Add lines of C_DMSP:C_tot?

# if (exportimg) {png(filename = paste0(opath,"dms_dmspt_",pg,"_counts.png"), width = 6, height = 7, units = 'cm', pointsize = 8, bg = 'white', res = 600, type = 'cairo')}
# 
# if (exportimg) {dev.off()}


# Plot for checking:
# plot(prof$Phaeo, 1e3*prof$phaeo_n_mL, xlab = "Phaeocystis cells/L, microscopy", ylab = "Phaeocystis cells/L, iFCB")
# plot(prof$Phaeo, 1e3*prof$phaeo_n_mL / prof$Phaeo, xlab = "Phaeocystis cells/L, microscopy", ylab = "Fraction detected: iFCB / microscopy")
# plot(prof$Phaeo_mgC_L.microscopy, prof$phaeo_mg_L, xlab = "Phaeocystis mgC/L, microscopy", ylab = "Phaeocystis mgC/L, iFCB")
plot(prof$diat_mgC_L.microscopy, prof$diat_pelagic_mg_L, xlab = "Diatoms mgC/L, microscopy", ylab = "Diatoms mgC/L, iFCB")
