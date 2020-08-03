# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(dplyr)
library(tidyr)
library(RColorBrewer)
library(classInt) # for function classIntervals
# Check http://geog.uoregon.edu/GeogR/topics/multiplots01.html
# detach("package:MASS", unload = T) # Interference with select function in dplyr

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
prof <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
diatL <- read.csv(file = "~/Desktop/GreenEdge/DiatomsLafond.tsv", header = T, sep = "\t")
diatL$diat_Laf.mgC_L.mic <- diatL$diat_Laf.mgC_L.mic/1000 # Correct units

# Exporting image?
exportimg <- F
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_phyto/"
ymax <- 158

# Use Cdmspt:Ctot ratio or intracellular concentration?
cmethod <- "ratio" # choose ratio or intraconc

addcd <- "_PhaeoCol_detr" # leave empty (default) or put _PhaeoCol_detr
if (cmethod == "intraconc"){
  addcd <- "" # Detritus and Phaeo colonies only allowed if using ratios
}

# Stefels2007 ratios? Does not apply if using quotas
Stefels <- "" # leave empty or put _Stefels
if(cmethod == "intraconc") {Stefels <- ""}

# Function to compute carbon content of phytoplankton
MendenDeuer_Lessard <- function(vol_microm3, is_diatom) {
  if (is_diatom) {
    pgC <- 0.288 * vol_microm3^0.811
  } else {
    pgC <- 0.216 * vol_microm3^0.939
  }
  return(pgC)
}

# ---------------------
# Data wrangling

# Merge Lafond diatoms data
prof <- merge(prof, diatL, by = c("stn","depth"))

# Correct units of biomass from ifcb
prof$detritus_mg_L <- prof$detritus_mg_L / 1000
prof$prym_clumped_mg_L <- prof$prym_clumped_mg_L / 1 # no correction needed
  
# Remove stations where microscopy counts not done
prof <- prof[,grep("NA",names(prof), invert = T)]
prof <- prof[prof$stn >= 418,]

# Remove data where no phyto counts available
prof <- prof[!is.na(prof$Phaeo) & !is.na(prof$dmspt) & !is.na(prof$depth),]
prof <- prof[prof$depth!=0.1 & !is.na(prof$prescreen) ,]

# Add surface and SCM categories
prof$scm <- 'SCM'
prof$scm[prof$depth < 10] <- 'surface'

# Rename DMS and select variables
cp <- prof$cpsmooth1
prof$dms <- prof$dms_consens_cf68
# prof <- select(prof, c(stn,depth,cast,year,month,day,scm,dmspt,diat_cen,diat_pen,dino_athec,dino_thec,chrys,crypt,Phaeo,flag,diat_Laf.mgC_L.mic,pras,prym_clumped_mg_L,detritus_mg_L))
prof <- select(prof, c(stn,depth,cast,year,month,day,scm,dmspt,diat_Laf.mgC_L.mic,dino_athec,dino_thec,Phaeo,prym_clumped_mg_L,chrys,crypt,flag,pras,detritus_mg_L,diat_cen,diat_pen))

# ---------------------
# Convert phyto counts to DMSP to estimate contribution

# Generalize to other species and make bar plots for surface and SCM. Calculate group-specific DMSPt and % in loop

# Assume DMSPp is 90% of DMSPt
fp <- 0.9

# Using Phaeocystis DMSP cell quota (q)
qdmsp.phaeo <- 10                                                 # fmol DMSP / cell
prof$dmspt.Phaeo.counts <- prof$Phaeo * qdmsp.phaeo * 1e-6        # DMSPt-Phaeo
prof$dmspt.Phaeo.counts.pc <- 100 * prof$dmspt.Phaeo.counts / (fp * prof$dmspt)
# print(summary(prof$dmspt.Phaeo.counts.pc))

# # Using Phaeocystis carbon from iFCB: WEIRD
# rdmsp.phaeo <- 0.05
# prof$dmspt.Phaeo.Cifcb <- prof$phaeo_mg_L * rdmsp.phaeo * (1/12) * (1/5) * 1e6
# prof$dmspt.Phaeo.Cifcb.pc <- 100*prof$dmspt.Phaeo.Cifcb/prof$dmspt
# print(summary(prof$dmspt.Phaeo.Cifcb.pc))

# Using biomass based on cell counts and ESD estimated from abundance-weighted average of size intervals (see Taxonomy xls v2018)
esd <- list(dino_athec = 15,
            dino_thec = 15,
            diat_cen = 20,
            diat_pen = 20,
            chrys = 10,
            crypt = 5,
            Phaeo = 5,
            flag = 5,
            pras = 14)
# # C_dmsp:C_tot ratio, optimized
rdmsp <- list(dino_athec = 0.11,
              dino_thec = 0.11,
              diat_cen = 0.02,
              diat_pen = 0.02,
              chrys = 0.10,
              crypt = 0.02,
              Phaeo = 0.1,
              flag = 0.05,
              diat_all = 0.02,
              pras = 0.025,
              detritus = 0.01)
if (Stefels == "_Stefels"){
# C_dmsp:C_tot ratio, Stefels2007
rdmsp <- list(dino_athec = 0.11,
              dino_thec = 0.11,
              diat_cen = 0.004,
              diat_pen = 0.004,
              chrys = 0.094,
              crypt = 0.02,
              Phaeo = 0.05,
              flag = 0.05,
              diat_all = 0.004,
              pras = 0.025,
              detritus = 0.01)}
# Intracellular DMSP conc millimolar. Now using rdmsp * 2000
cdmsp <- list(dino_athec = 200, 
              dino_thec = 200,
              diat_cen = 20,
              diat_pen = 20,
              chrys = 200,
              crypt = 40,
              Phaeo = 100,
              flag = 100,
              diat_all = 20,
              pras = 50)

for (nn in names(esd)) {
  
  vol_microm3 <- (4/3) * pi * ( (esd[[nn]] /2)^3 ) # µm^3
  
  if (cmethod == "ratio") {
    
    # Option A: from cell carbon and C_dmsp:C_tot ratio, Menden-Deuer and Lessard
    if (nn %in% c("diat_cen","diat_pen")) {
      tmp <- MendenDeuer_Lessard(vol_microm3, is_diatom = T) * 1e-9
    } else {
      tmp <- MendenDeuer_Lessard(vol_microm3, is_diatom = F) * 1e-9
    }
    prof[[paste(nn,"mgC_L.mic",sep = ".")]] <- tmp * prof[[nn]]
    prof[[paste("dmspt",nn,"Cmic",sep = ".")]] <- tmp * prof[[nn]] * rdmsp[[nn]] * (1/12) * (1/5) * 1e6
    prof[[paste("dmspt",nn,"Cmic.pc",sep = ".")]] <- 100 * tmp * prof[[nn]] * rdmsp[[nn]] * (1/12) * (1/5) * 1e6 / (fp * prof$dmspt)
    
  } else if (cmethod == "intraconc") {
    
    # Option B: from cell volume
    prof[[paste("dmspt",nn,"Cmic",sep = ".")]] <- vol_microm3 * prof[[nn]] * cdmsp[[nn]] * (1e6 / 1e15)
    prof[[paste("dmspt",nn,"Cmic.pc",sep = ".")]] <- 100 * vol_microm3 * prof[[nn]] * cdmsp[[nn]] * (1e6 / 1e15) / (fp * prof$dmspt)
  }
}

if (cmethod == "ratio") {
  # # Compare Lafond estimates and ESD-based estimates. ESD of 20 for both centric and pennates gives good approximation
  # plot(prof$diat_Laf.mgC_L.mic, prof$diat_cen.mgC_L.mic + prof$diat_pen.mgC_L.mic,
  #      main = paste0("r = ",cor(prof$diat_Laf.mgC_L.mic, prof$diat_cen.mgC_L.mic + prof$diat_pen.mgC_L.mic, "pairwise.complete.obs", method = "s")))
}

# Remove microscopy counts-based data for diatoms
if (cmethod == "ratio") {
  # prof <- select(prof, -c(diat_cen.mgC_L.mic, dmspt.diat_cen.Cmic, dmspt.diat_cen.Cmic.pc, diat_pen.mgC_L.mic, dmspt.diat_pen.Cmic, dmspt.diat_pen.Cmic.pc,dmspt.Phaeo.counts,dmspt.Phaeo.counts.pc))
  prof <- select(prof, -c(diat_cen.mgC_L.mic, dmspt.diat_cen.Cmic, dmspt.diat_cen.Cmic.pc, diat_pen.mgC_L.mic, dmspt.diat_pen.Cmic, dmspt.diat_pen.Cmic.pc))  
} else if (cmethod == "intraconc") {
  prof <- select(prof, -c( dmspt.diat_cen.Cmic, dmspt.diat_cen.Cmic.pc, dmspt.diat_pen.Cmic, dmspt.diat_pen.Cmic.pc))
}

# Add DMSP for Lafond diatoms
prof$dmspt.diat_Laf.Cmic <- prof$diat_Laf.mgC_L.mic * rdmsp$diat_all * (1/12) * (1/5) * 1e6
prof$dmspt.diat_Laf.Cmic.pc <- 100 * prof$dmspt.diat_Laf.Cmic / (fp * prof$dmspt)

if (addcd == "_PhaeoCol_detr"){
  # Add DMSP for colonies "clumped prym". Assume 1% of colonies contain DMSP and rdmsp = 0.10 (equivalent: 2% and 0.05)
  prof$dmspt.Phaeo_col.ifcb <- (prof$prym_clumped_mg_L / 100) * rdmsp$Phaeo * (1/12) * (1/5) * 1e6
  prof$dmspt.Phaeo_col.ifcb.pc <- 100 * prof$dmspt.Phaeo_col.ifcb / (fp * prof$dmspt)
  # Add DMSP for detritus
  prof$dmspt.detritus.ifcb <- prof$detritus_mg_L * rdmsp$detritus * (1/12) * (1/5) * 1e6 # NOTE: units corrected from mg to ugC/L (above)
  prof$dmspt.detritus.ifcb.pc <- 100 * prof$dmspt.detritus.ifcb / (fp * prof$dmspt)
  # Increase maximum y axis DMSP
  ymax <- 200
}

# OUTPUT
OUT <- prof
# View(OUT)
write.csv(x = OUT, file = paste0(opath,"Fraction_DMSPt_phyto",addcd,".csv"), row.names = F)


# ---------------------------------------------------------------
# Plot

# Bar color scheme
# if (addcd == "_PhaeoCol_detr"){
  phycol <- brewer.pal(12, 'Paired')
# } else {
  # phycol <- brewer.pal(8, 'Dark2')[c(1,2,3,4,5,7,8,6)]
# }

phynames <- list("dmspt.diat_Laf.Cmic.pc"="Diat_all",
                 "dmspt.diat_cen.Cmic.pc"="Diat_C",
                 "dmspt.diat_pen.Cmic.pc"="Diat_P",
                 "dmspt.dino_athec.Cmic.pc"="Dino_A",
                 "dmspt.dino_thec.Cmic.pc"="Dino_T",
                 "dmspt.chrys.Cmic.pc"="Chryso",
                 "dmspt.crypt.Cmic.pc"="Crypto",
                 "dmspt.Phaeo.Cmic.pc"="Phaeocystis",
                 "dmspt.Phaeo.counts.pc"="Phaeocystis_cq",
                 "dmspt.flag.Cmic.pc"="Flag_other",
                 "dmspt.Phaeo_col.ifcb.pc"="Phaeo_col",
                 "dmspt.detritus.ifcb.pc"="Detritus",
                 "dmspt.pras.Cmic.pc"="Prasino")
dfplot <-prof[ , c( "stn","depth","scm","dmspt",names(prof)[grepl("dmspt.", names(prof)) & grepl(".pc", names(prof))] )]
for (nn in names(phynames)) {names(dfplot)[names(dfplot) == nn] <- phynames[[nn]]}
for (nn in names(phynames)) {print(phynames[[nn]])}

# # View sum of percentages
# sums <- data.frame('stn' = dfplot[dfplot$scm=="surface", "stn"],
#                    'surface' = rowSums(dfplot[dfplot$scm=="surface", 6:dim(dfplot)[2]]),
#                    'SCM' = rowSums(dfplot[dfplot$scm=="SCM", 6:dim(dfplot)[2]]))
# View(sums)
# # TidyR format: not needed if using base barplot
# gather(prof[ ,grepl("dmspt.", names(prof)) & !grepl(".pc", names(prof))])


if (exportimg) {png(filename = paste0(opath,"dmspt_biomass_v3_",cmethod,addcd,Stefels,".png"), width = 8, height = 12, units = 'cm', pointsize = 8, bg = 'white', res = 600, type = 'cairo')}

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
        col = phycol,
        border = F,
        ylab = "% of DMSPp",
        # ylab = "",
        axes = F,
        ylim = c(0,ymax),
        cex = 1.2,
        cex.names = 1.2,
        cex.lab = 1.2,
        las = 1,
        main = "a) Surface")
        # main = "b) Surface, adjusted")
abline(h = 100, col = "gray80", lwd = 1)
mtext(text = as.character(round(fp*dfplot[dfplot$scm=="surface", "dmspt"], digits = 0)),
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
        col = phycol,
        border = F,
        ylab = "% of DMSPp",
        # ylab = "",
        axes = F,
        ylim = c(0,ymax),
        cex = 1.2,
        cex.names = 1.2,
        cex.lab = 1.2,
        las = 1,
        main = "b) SCM",
        # main = "d) SCM, adjusted",
        xlab = "Station")
abline(h = 100, col = "gray80", lwd = 1)
mtext(text = as.character(round(fp*dfplot[dfplot$scm=="SCM", "dmspt"], digits = 0)),
      side = 3,
      line = -3,
      at = seq(0.6, 0.6+8*1.2, 1.2),
      cex = 0.8,
      font = 2)
axis(2, labels = T, cex.axis = 1.2)

if (exportimg) {dev.off()}


# ---------------------------------------------------------------
# Dino and diatom DMSPp % contribution
print(sapply(dfplot[,c("Dino_A","Diat_all")], mean, na.rm=T))
print(sapply(dfplot[,c("Dino_A","Diat_all")], sd, na.rm=T))

# ---------------------------------------------------------------
# Agreement between indirect estimates and measurements

if (addcd == "_PhaeoCol_detr"){
  dmsppSUM <- rowSums(dfplot[ , c("Dino_A","Dino_T","Chryso","Crypto","Phaeocystis","Flag_other","Prasino","Diat_all","Phaeo_col","Detritus") ], na.rm = T)
} else {
  dmsppSUM <- rowSums(dfplot[ , c("Dino_A","Dino_T","Chryso","Crypto","Phaeocystis","Flag_other","Prasino","Diat_all") ], na.rm = T)  
}

print( mean(dmsppSUM[dfplot$scm=="surface"]) )
print( range(dmsppSUM[dfplot$scm=="surface"]) )

print( mean(dmsppSUM[dfplot$scm=="SCM"]) )
print( range(dmsppSUM[dfplot$scm=="SCM"]) )

# ---------------------------------------------------------------
# Checks

# # Check relationship between Phaeocystis solitary and colonies, by SCM vs. surface
# plot(prof$Phaeo.mgC_L.mic, prof$prym_clumped_mg_L)
# points(prof$Phaeo.mgC_L.mic[prof$scm=='SCM'], prof$prym_clumped_mg_L[prof$scm=='SCM'], pch = 19)
# 
# plot(prof$Phaeo.mgC_L.mic, prof$prym_clumped_mg_L/prof$Phaeo.mgC_L.mic)


# Check relationship between POC estimated from Cp and total carbon biomass from my estimates, Lafond and iFCB
cf <- 391 # [mg C/m2], Cetinic et al. 2012
prof$poc <- cp * cf # [mg C/m3] 

# plot(prof$poc, prof$diat_Laf.mgC_L.mic)



