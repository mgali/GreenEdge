# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(dplyr)
library(RColorBrewer)
library(classInt) # for function classIntervals
# Check http://geog.uoregon.edu/GeogR/topics/multiplots01.html

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
prof <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
diatL <- read.csv(file = "~/Desktop/GreenEdge/DiatomsLafond.tsv", header = T, sep = "\t")
diatL$diat_Laf.mgC_L.mic <- diatL$diat_Laf.mgC_L.mic/1000 # Correct units

# Exporting image?
exportimg <- T
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_phyto/"

# Function to compute carbon content of phytoplankton
MendenDeuer_Lessard <- function(esd_microm, cell_conc, is_diatom) {
  vol <- (4/3) * pi * ((esd_microm/2)^3) # µm^3
  if (is_diatom) {
    pgC <- 0.288 * vol^0.811
  } else {
    pgC <- 0.216 * vol^0.939
  }
  return(pgC * cell_conc)
}

# ---------------------
# Data wrangling

# Merge Lafond diatoms data
prof <- merge(prof, diatL, by = c("stn","depth"))

# Remove stations where microscopy counts not done
prof <- prof[,grep("NA",names(prof), invert = T)]
prof <- prof[prof$stn >= 418,]

# Remove data where no phyto counts available
prof <- prof[!is.na(prof$Phaeo) & !is.na(prof$dmspt) & !is.na(prof$depth),]
prof <- prof[prof$depth!=0.1 & !is.na(prof$prescreen) ,]

# Add surface and SCM categories
prof$scm <- 'SCM'
prof$scm[prof$depth <= 10] <- 'surface'

# Rename DMS and remove unnecessary variables
prof$dms <- prof$dms_consens_cf68
prof <- select(prof, c(stn,depth,cast,year,month,day,scm,dmspt,diat_cen,diat_pen,dino_athec,dino_thec,chrys,crypt,Phaeo,flag,diat_Laf.mgC_L.mic))

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
esd <- list(dino_athec = 15,
            dino_thec = 15,
            diat_cen = 20,
            diat_pen = 20,
            chrys = 10,
            crypt = 5,
            Phaeo = 5,
            flag = 5)
rdmsp <- list(dino_athec = 0.10,
              dino_thec = 0.10,
              diat_cen = 0.02,
              diat_pen = 0.02,
              chrys = 0.10,
              crypt = 0.02,
              Phaeo = 0.10,
              flag = 0.10,
              diat_all = 0.02)

for (nn in names(esd)) {
  if (nn %in% c("diat_cen","diat_pen")) {
    tmp <- MendenDeuer_Lessard(esd_microm = esd[[nn]], cell_conc = prof[[nn]], is_diatom = T) * 1e-9
  } else {
    tmp <- MendenDeuer_Lessard(esd_microm = esd[[nn]], cell_conc = prof[[nn]], is_diatom = F) * 1e-9
  }
  prof[[paste(nn,"mgC_L.mic",sep = ".")]] <- tmp
  prof[[paste("dmspt",nn,"Cmic",sep = ".")]] <- tmp * rdmsp[[nn]] * (1/12) * (1/5) * 1e6
  prof[[paste("dmspt",nn,"Cmic.pc",sep = ".")]] <- 100 * tmp * rdmsp[[nn]] * (1/12) * (1/5) * 1e6 / (fp * prof$dmspt)
}

# # Compare Lafond estimates and ESD-based estimates. ESD of 20 for both centric and pennates gives good approximation
# plot(prof$diat_Laf.mgC_L.mic, prof$diat_cen.mgC_L.mic + prof$diat_pen.mgC_L.mic,
#      main = paste0("r = ",cor(prof$diat_Laf.mgC_L.mic, prof$diat_cen.mgC_L.mic + prof$diat_pen.mgC_L.mic, "pairwise.complete.obs", method = "s")))

# Add DMSP for Lafond diatoms
prof$dmspt.diat_Laf.Cmic <- prof$diat_Laf.mgC_L.mic * rdmsp$diat_all * (1/12) * (1/5) * 1e6
prof$dmspt.diat_Laf.Cmic.pc <- 100 * prof$dmspt.diat_Laf.Cmic / (fp * prof$dmspt)

# Remove microscopy counts-based data
prof <- select(prof, -c(diat_cen.mgC_L.mic, dmspt.diat_cen.Cmic, dmspt.diat_cen.Cmic.pc, diat_pen.mgC_L.mic, dmspt.diat_pen.Cmic, dmspt.diat_pen.Cmic.pc))

# Add sum of percentages
prof$sum.pc <- rowSums(prof[ , grepl(".pc",names(prof))])
# View(prof$sum.pc)

# OUTPUT
OUT <- prof
# View(OUT)
write.csv(x = OUT, file = paste0(opath,"Fraction_DMSPt_phyto.csv"), row.names = F)


# ---------------------
# Plot

phynames <- list("dmspt.diat_Laf.Cmic"="Diat_all",
                 "dmspt.diat_cen.Cmic"="Diat_C",
                 "dmspt.diat_pen.Cmic"="Diat_P",
                 "dmspt.dino_athec.Cmic"="Dino_A",
                 "dmspt.dino_thec.Cmic"="Dino_T",
                 "dmspt.chrys.Cmic"="Chryso",
                 "dmspt.crypt.Cmic"="Crypto",
                 "dmspt.Phaeo.Cmic"="Phaeocystis",
                 "dmspt.flag.Cmic"="Flag_other")
# dfplot <- 


# if (exportimg) {png(filename = paste0(opath,"dms_dmspt_",pg,"_counts.png"), width = 6, height = 7, units = 'cm', pointsize = 8, bg = 'white', res = 600, type = 'cairo')}

# Multipanel setup
# m <- rbind(matrix(data = 1, nrow = 4, ncol = 9), matrix(data = 2, nrow = 4, ncol = 9))
# layout(m)
# par(oma = c(1,1,0.5,0.5))

# ---------------------
# a) DMSPp partition surface


# ---------------------
# b) DMSPp partition SCM


# if (exportimg) {dev.off()}



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