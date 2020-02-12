# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(dplyr)
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

# Remove stations where microscopy counts not done
prof <- prof[,grep("NA",names(prof), invert = T)]
prof <- prof[pplot$stn >= 418,]

# Remove data where no phyto counts available
prof <- prof[!is.na(prof$Phaeo) & !is.na(prof$dmspt) & !is.na(prof$depth),]
prof <- prof[prof$depth!=0.1 & !is.na(prof$prescreen) ,]

# Add surface and SCM categories
prof$scm <- 'SCM'
prof$scm[prof$depth <= 10] <- 'surface'

# Rename DMS and remove unnecessary variables
prof$dms <- prof$dms_consens_cf68
prof <- select(prof, c(stn,depth,cast,year,month,day,scm,dmspt,diat_cen,diat_pen,dino_athec,dino_thec,chrys,crypt,Phaeo,flag))

# ---------------------
# Convert phyto counts to DMSP to estimate contribution

# Generalize to other species and make bar plots for surface and SCM. Calculate group-specific DMSPt and % in loop

# Assume DMSPp is 90% of DMSPt
fp <- 0.9

# Using Phaeocystis DMSP cell quota (q)
qdmsp.phaeo <- 10                                                 # fmol DMSP / cell
prof$dmspt.Phaeo.counts <- prof$Phaeo * qdmsp.phaeo * 1e-6        # DMSPt-Phaeo
prof$dmspt.Phaeo.counts.pc <- 100 * prof$dmspt.Phaeo.counts / (fp * prof$dmspt)
print(summary(prof$dmspt.Phaeo.counts.pc))

# # Using Phaeocystis carbon from iFCB: WEIRD
# rdmsp.phaeo <- 0.05
# prof$dmspt.Phaeo.Cifcb <- prof$phaeo_mg_L * rdmsp.phaeo * (1/12) * (1/5) * 1e6
# prof$dmspt.Phaeo.Cifcb.pc <- 100*prof$dmspt.Phaeo.Cifcb/prof$dmspt
# print(summary(prof$dmspt.Phaeo.Cifcb.pc))

# Using carbon biomass based on cell counts and ESD estimated from abundance-weighted average of size intervals (see Taxonomy xls v2018)
rdmsp <- list(diat_cen = 0.01,
              diat_pen = 0.03,
              dino_athec = 0.10,
              dino_thec = 0.10,
              chrys = 0.10,
              crypt = 0.02,
              Phaeo = 0.10,
              flag = 0.10)
esd <- list(diat_cen = 100,
            diat_pen = 20,
            dino_athec = 15,
            dino_thec = 15,
            chrys = 10,
            crypt = 5,
            Phaeo = 5,
            flag = 5)

for (nn in names(rdmsp)) {
  if (nn %in% c("diat_cen","diat_pen")) {
    tmp <- MendenDeuer_Lessard(esd_microm = esd[[nn]], cell_conc = prof[[nn]], is_diatom = T) * 1e-9
  } else {
    tmp <- MendenDeuer_Lessard(esd_microm = esd[[nn]], cell_conc = prof[[nn]], is_diatom = F) * 1e-9
  }
  prof[[paste(nn,"mgC_L.mic",sep = ".")]] <- tmp
  prof[[paste("dmspt",nn,"Cmic",sep = ".")]] <- tmp * rdmsp[[nn]] * (1/12) * (1/5) * 1e6
  prof[[paste("dmspt",nn,"Cmic.pc",sep = ".")]] <- 100 * tmp * rdmsp[[nn]] * (1/12) * (1/5) * 1e6 / (fp * prof$dmspt)
}

# OUTPUT
OUT <- prof
View(OUT)
write.csv(x = OUT, file = paste0(opath,"Fraction_DMSPt_phyto.csv"), row.names = F)


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


# --------------------------
# OLD
# prof$Phaeo_mgC_L.microscopy <- MendenDeuer_Lessard(esd_microm = 5, cell_conc = prof$Phaeo, is_diatom = F) * 1e-9
# rdmsp.phaeo <- 0.10
# prof$dmspt_phaeo.Cmic <- prof$Phaeo_mgC_L.microscopy * rdmsp.phaeo * (1/12) * (1/5) * 1e6
# prof$dmspt_phaeo.Cmic.pc <- 100 * prof$dmspt_phaeo.Cmic / (fp * prof$dmspt)
# print(summary(prof$dmspt_phaeo.Cmic.pc))
# 
# prof$diat_cen_mgC_L.microscopy <- MendenDeuer_Lessard(esd_microm = 100, cell_conc = prof$diat_cen, is_diatom = T) * 1e-9
# prof$diat_pen_mgC_L.microscopy <- MendenDeuer_Lessard(esd_microm = 20, cell_conc = prof$diat_pen, is_diatom = T) * 1e-9
# prof$diat_mgC_L.microscopy <- prof$diat_cen_mgC_L.microscopy + prof$diat_pen_mgC_L.microscopy
# rdmsp.diat <- 0.02
# prof$dmspt_diat.Cmic <- prof$diat_mgC_L.microscopy * rdmsp.diat * (1/12) * (1/5) * 1e6
# prof$dmspt_diat.Cmic.pc <- 100 * prof$dmspt_diat.Cmic / (fp * prof$dmspt)
# print(summary(prof$dmspt_diat.Cmic.pc))

# OUT <- prof[,c("stn","depth",
#                "Phaeo","diat_cen","diat_pen",
#                "Phaeo_mgC_L.microscopy","diat_cen_mgC_L.microscopy","diat_pen_mgC_L.microscopy","diat_mgC_L.microscopy",
#                "dmspt.Phaeo.counts","dmspt_phaeo.Cmic","dmspt_diat.Cmic",
#                "dmspt",
#                "dmspt.Phaeo.counts.pc","dmspt_phaeo.Cmic.pc","dmspt_diat.Cmic.pc")]