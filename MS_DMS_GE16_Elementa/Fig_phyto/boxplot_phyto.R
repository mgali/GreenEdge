# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(dplyr)
library(ggplot2)
library(reshape2) # to convert wide to long using melt

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
surf.all <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)
FC <- read.csv(file = "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_phyto/flow_cyto_preprocessed.Rda", header = T) # from flo_cyto_preprocess.R

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

# Remove rows where phyto counts are missing
pplot <- pplot[!is.na(pplot$Phaeo) & !is.na(pplot$diat_cen),]

# Remove duplicated rows
dd <- duplicated(pplot[,c("cast","depth")]) | (pplot$station==418 & pplot$depth==0.1)
pplot <- pplot[!dd,]

# Add surface and SCM categories
pplot$scm <- 'SCM'
pplot$scm[pplot$depth < 10] <- 'surface'
# pplot$scm <- 'Subsurface Chl max'
# pplot$scm[pplot$depth < 10] <- 'Surface'


# Merge flow cytometry
pplot <- merge(x = pplot, y = FC[,c("cast","scm","Pico_FC","Nano_FC","Crypto_FC","Bact_FC")], by = c("cast","scm"), all.x = T)

xl <- list("diat_cen"="Diatoms (centric)",
           "diat_pen"="Diatoms (pennate)",
           "dino_athec"="Dinoflagellates (athecate)",
           "dino_thec"="Dinoflagellates (thecate)",
           # "chlor"="Chlorophytes",
           "dictyo"="Dictyochophytes",
           "chrys"="Chrysophytes",
           "crypt"="Cryptophytes",
           "pras"="Prasino. (Pyramimonas)",
           "Pico_FC"="Picophytoeukaryotes (other)",
           "flag"="Nanoflagellates (other)",
           # "prym"="Prymnesiophyceae", # Nearly equivalent to Phaeocystis
           # "Phaeo"=expression(italic("Phaeocystis pouchetii"))) # italics with expression() not working here
           "Phaeo"="Phaeocystis pouchetii")
pplot[,names(xl)] <- pplot[,names(xl)]+1

pplot <- melt(data = pplot[,c("station","scm",names(xl))],
              measure.vars = names(xl),
              value.name = "counts") # Abundance (cells/mL)
pplot$group <- factor(pplot$variable)
pplot$Depth <- factor(pplot$scm)
# pplot$Depth <- factor(pplot$Depth, levels = rev(levels(pplot$Depth)))

if (exportimg) {png(filename = paste0(opath,"FigS1_phytoCounts.png"), width = 16, height = 8, units = 'cm', pointsize = 8, bg = 'white', res = 600, type = 'cairo')}

# Phyto counts NEED TO ADD PEUK
p <- ggplot(pplot, aes(x=group, y=counts, fill=Depth)) +
  scale_y_continuous(trans='log10', limits = c(2e2,2e7)) +
  # scale_x_discrete(limits = rev(levels(pplot$Depth))) +
  scale_x_discrete(labels = xl) +
  coord_flip() +
  geom_boxplot(notch = F,
               outlier.shape = 19,
               outlier.size = 0.4,
               outlier.stroke = 0.5) + 
  # scale_fill_brewer(palette="YlGnBu") +
  scale_fill_manual(values = rev(brewer.pal(4,"YlGnBu"))[c(2,4)]) +
  # geom_jitter(aes(colour = Depth), width = 0.1, size = 0.1) +
  xlab("") +
  ylab(expression(paste("Abundance ","(cells ",L^-1*")"))) +
  theme_bw()
print(p)

if (exportimg) {dev.off()}

