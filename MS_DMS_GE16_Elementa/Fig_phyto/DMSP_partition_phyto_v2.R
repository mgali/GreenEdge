# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(dplyr)
library(tidyr)
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
prof$scm[prof$depth < 10] <- 'surface'

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
rdmsp <- list(dino_athec = 0.05,
              dino_thec = 0.05,
              diat_cen = 0.01,
              diat_pen = 0.01,
              chrys = 0.10,
              crypt = 0.02,
              Phaeo = 0.10,
              flag = 0.10,
              diat_all = 0.01)

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

# Remove microscopy counts-based data for diatoms
prof <- select(prof, -c(diat_cen.mgC_L.mic, dmspt.diat_cen.Cmic, dmspt.diat_cen.Cmic.pc, diat_pen.mgC_L.mic, dmspt.diat_pen.Cmic, dmspt.diat_pen.Cmic.pc))

# OUTPUT
OUT <- prof
# View(OUT)
write.csv(x = OUT, file = paste0(opath,"Fraction_DMSPt_phyto.csv"), row.names = F)

# ---------------------
# Plot

parcol <- brewer.pal(8, 'Dark2')

phynames <- list("dmspt.diat_Laf.Cmic.pc"="Diat_all",
                 "dmspt.diat_cen.Cmic.pc"="Diat_C",
                 "dmspt.diat_pen.Cmic.pc"="Diat_P",
                 "dmspt.dino_athec.Cmic.pc"="Dino_A",
                 "dmspt.dino_thec.Cmic.pc"="Dino_T",
                 "dmspt.chrys.Cmic.pc"="Chryso",
                 "dmspt.crypt.Cmic.pc"="Crypto",
                 "dmspt.Phaeo.Cmic.pc"="Phaeocystis",
                 "dmspt.Phaeo.counts.pc"="Phaeocystis_cq",
                 "dmspt.flag.Cmic.pc"="Flag_other")
dfplot <-prof[ , c( "stn","depth","scm","dmspt",names(prof)[grepl("dmspt.", names(prof)) & grepl(".pc", names(prof))] )]
for (nn in names(phynames)) {names(dfplot)[names(dfplot) == nn] <- phynames[[nn]]}

# # View sum of percentages
# sums <- data.frame('stn' = dfplot[dfplot$scm=="surface", "stn"],
#                    'surface' = rowSums(dfplot[dfplot$scm=="surface", 6:dim(dfplot)[2]]),
#                    'SCM' = rowSums(dfplot[dfplot$scm=="SCM", 6:dim(dfplot)[2]]))
# View(sums)

# # TidyR format
# gather(prof[ ,grepl("dmspt.", names(prof)) & !grepl(".pc", names(prof))])


if (exportimg) {png(filename = paste0(opath,"dmspt_biomass.png"), width = 16, height = 12, units = 'cm', pointsize = 8, bg = 'white', res = 600, type = 'cairo')}

# Multipanel setup
m <- rbind(matrix(data = 1, nrow = 4, ncol = 9), matrix(data = 2, nrow = 4, ncol = 9))
layout(m)
par(oma = c(1,1,0.5,0.5))

# ---------------------
# a) DMSPp partition surface
par(mar = c(5,5,2,2))

barplot(height = t(as.matrix(dfplot[dfplot$scm=="surface",6:dim(dfplot)[2]])),
        names.arg = rep("",9),
        beside = F,
        legend.text = names(dfplot)[6:dim(dfplot)[2]],
        args.legend = list(x = "right", bty = "n", cex = 1.2),
        col = parcol,
        border = F,
        ylab = "% of DMSPp",
        axes = F,
        ylim = c(0,100),
        cex = 1.2,
        cex.names = 1.2,
        cex.lab = 1.2,
        las = 1,
        main = "a) Surface")
mtext(text = as.character(round(dfplot[dfplot$scm=="surface", "dmspt"], digits = 0)),
      side = 3,
      line = -3,
      at = seq(0.6, 0.6+8*1.2, 1.2),
      cex = 0.8,
      font = 2)
axis(2, labels = T, cex.axis = 1.2)

# ---------------------
# b) DMSPp partition SCM
par(mar = c(5,5,2,2))

barplot(height = t(as.matrix(dfplot[dfplot$scm=="SCM",6:dim(dfplot)[2]])),
        names.arg = dfplot$stn[dfplot$scm=="SCM"],
        beside = F,
        col = parcol,
        border = F,
        ylab = "% of DMSPp",
        axes = F,
        ylim = c(0,100),
        cex = 1.2,
        cex.names = 1.2,
        cex.lab = 1.2,
        las = 1,
        main = "b) Subsurface chlorophyll maximum (SCM)",
        xlab = "Station")
mtext(text = as.character(round(dfplot[dfplot$scm=="SCM", "dmspt"], digits = 0)),
      side = 3,
      line = -3,
      at = seq(0.6, 0.6+8*1.2, 1.2),
      cex = 0.8,
      font = 2)
axis(2, labels = T, cex.axis = 1.2)

if (exportimg) {dev.off()}

