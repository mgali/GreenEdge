# ESTIMATE BACYERIAL DMSP STOCK

library("readxl")
library(RColorBrewer)
library(dplyr)

# Load data
sulfur.all <- read.csv(file = "~/Desktop/GreenEdge/GCMS/GE2016.profiles.ALL.OK.csv")
prof.all <- read.csv(file = "~/Desktop/GreenEdge/Achim_Randelhoff/data/Randelhoff-et-al_GreenEdge_ANP_v0.1.csv", header = T)
surf.all <- read.csv(file = "~/Desktop/GreenEdge/Achim_Randelhoff/data/Randelhoff-et-al_Ice-edge-paper_per-station_v0.1.csv", header = T)
pre.qpcr <- read_excel(path = "~/Desktop/GreenEdge/Results_methylotrophy story_ODV2-1_mg.xlsx", sheet = "Results_methylotrophy story_ODV")

# Exporting image?
# exportimg <- T # Figure with binned profiles by OWD categories
# opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_qPCR_norm_grouped/"

# ---------------------
# Merge with profiles to get clustering coefficient and SIC classification
pplot <- merge(x = prof.all, y = surf.all, all.x = T, all.y = F, by = 'Station', suffixes = "")

# Preprocess qpcr data
qpcr <- data.frame(Station = pre.qpcr$Station,
                   cast = pre.qpcr$CTD,
                   Depth_m = pre.qpcr$Depth,
                   dddP = pre.qpcr$`Average of dddP copies/mL`,
                   dmdA = pre.qpcr$`Average of dmdA copies /mL`,
                   Methylobacterium = pre.qpcr$Methylobacterium*100,
                   Rhodobacteraceae = pre.qpcr$Rhodobacteraceae*100,
                   SAR11  = pre.qpcr$SAR11*100,
                   SAR116 = pre.qpcr$SAR116*100,
                   Polaribacter  = pre.qpcr$Polaribacter*100,
                   Thiotrichales  = pre.qpcr$Thiotrichales*100,
                   Methylophylales = pre.qpcr$Methylophylales*100,
                   Oceanospirillales = pre.qpcr$Oceanospirillales*100,
                   Pseudoalteromonas = pre.qpcr$Pseudoalteromonas*100,
                   BA = pre.qpcr$`BA cells/mL`,
                   BP = pre.qpcr$`Bacterial production µgC/L/d`
)
qpcr$Station <- as.numeric(gsub(pattern = "G", replacement = "", qpcr$Station))

# Merge with qpcr
pplot <- merge(x = qpcr, y = pplot, all.x = T, all.y = F, by = c('Station','Depth_m'), suffixes = "")

# Merge with main sulfur profiles data
psulfur <- data.frame(dmspt = sulfur.all$dmspt,
                      dms = sulfur.all$dms_consens_cf68,
                      Station = sulfur.all$stn,
                      Depth_m = sulfur.all$depth)
pplot <- merge(x = psulfur, y = pplot, all.x = T, all.y = F, by = c('Station','Depth_m'), suffixes = "")

# Remove duplicated columns with NA in their names
pplot <- pplot[,grep("NA",names(pplot), invert = T)]

# Remove duplicated rows 
dd <- duplicated(pplot[,c("cast","Depth_m")])
pplot <- pplot[!dd,]

# # Hide data from stn 400 (either entire or just surface)
# pplot[pplot$Station>=400,] <- NA

# ------------------------------------------------------------------------------------
# Estimate bacterial biovolume --> protein --> DMSP
bp <- 4.8 * 1e-9                                # bact protein (µg == fg 1e-9) for a bacterium of 0.51 µm diameter (Simon & Azam, 1989)
dmsp2bp <- 40                                   # pmol DMSP/µg protein, N starved bacterium (upper bound; Curson et al., 2017)
pplot$Bdmsp <- pplot$BA * bp * dmsp2bp          # pmol_DMSP/mL == nM; cells/mL * µg_prot/cell * dmsp/µg_prot


# Plots, printouts
plot(pplot$dmspt, pplot$Bdmsp, type="p")
print(summary(pplot$Bdmsp))
print(summary(100*pplot$Bdmsp/pplot$dmspt))
