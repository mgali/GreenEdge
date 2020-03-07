# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(dplyr)
library(MASS)
library(corrplot) # install.packages("corrplot")
# library(Hmisc) #install.packages("Hmisc")

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
prof <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)

# Output settings
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/regression_stepwise/"
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

# Count proportion of empty cells in each column (mostly used to exclude some pigments or group them into functional units as DD or VAZ xanthophylls)
noNAcount <- sapply(prof, function(x) round( sum(!is.na(x) & (!is.na(prof$dms) | !is.na(prof$dmspt))), digits = 2) )

# Group DD and VAZ cycles, PSC and PPC
prof$dd <- rowSums(prof[,c("diadino","diato")], na.rm = T)
prof$vaz <- rowSums(prof[,c("zea","anthera","viola")], na.rm = T)
prof$psc <- rowSums(prof[,c("fuco","peri","but19_like","hex","hex19_likeSUM")], na.rm = T)
prof$ppc <- rowSums(prof[,c("zea","anthera","viola","diadino","diato","allo","tcar")], na.rm = T)

# ------------------------------------------------------------------------
# Prepare regression and run in loop
# Pigments with less than 50% available data

yvarS <- c("dmspt","dms")
xvarS <- list(physics = c("temp","sal","N2","cpsmooth1","anp","par_d_p24h_ein_m_2_day_1","dmspt"),
              pigments = c("chlc3","chlc2group","chldaSUM","peri","phdaSUM","but","fuco","neo","pras.1","hex",
                           "dd","allo","lut","chlb","tchla","phytnaSUM","tcar","but19_like"),
              # pigments_reduced = c("chlc3","chlc2group","peri","neo","hex","dd","allo","chlb","tchla","tcar","but19_like"),
              # pigments_reduced = c("chlc3","tchla","but19_like"), # Sample size is 56 due to missing values
              # pigments_reduced = c("chlc3","chlc2group","psc","ppc"),
              # pigments4dmsp = c("chlc2group","hex","tchla","tcar","fuco")) # Best subset for DMSPt
              pigments4dmspBIS = c("but19_like","chlc2group","hex","tchla","fuco")) # Best subset for DMSPt BIS: idem adding 19'but-like. and removing tcar
              # pigments4dms = c("but19_like","chlc2group","hex","chlb","tcar","fuco")) # Best subset for DMS. Replacing tchla by chlb makes a difference. Adding peri: marginally significant and improves R2 to 0.82
              # pigments4dmsBIS = c("but19_like","chlc2group","hex","chlb","tcar","fuco")) # Best subset for DMS BIS: idem adding 19'but-like. Very high R2 but N drops by 1/3

for (mm in c("spearman","pearson")) {
  for (yvar in yvarS[2]) { #[1]
    for (xi in names(xvarS)[3]) { #[2]
      
      profOK <- prof
      if (xi == "pigments") { profOK <- prof[prof$stn != 413, ] }
      XY <- profOK[, c(yvar,xvarS[[xi]])]
      if (xi == "pigments") {
        TMP <- XY[,2:dim(XY)[2]]
        # TMP[is.na(TMP)] <- 0
        XY[,2:dim(XY)[2]] <- TMP
      }
      
      # # Correlation matrix: save plot and data
      # oname <- paste(yvar,xi,mm, sep = "_")
      # res <- round(cor(XY, method = mm, use = "pairwise"), 2)
      # rdf <- data.frame(yvar = res[ , colnames(res) == yvar])
      # # View(rdf)
      # if (exportimg) {
      #   png(filename = paste0(opath,oname,".png"), width = 17, height = 17, units = 'cm', pointsize = 8, bg = 'white', res = 600, type = 'cairo')
      # }
      # corrplot(res, type = "upper", order = "hclust",
      #          tl.col = "black", tl.srt = 45)
      # if (exportimg) {dev.off()}
      # write.csv(rdf, file = paste0(opath,oname,".csv"), row.names = T)
      
      # Stepwise regression
      XY <- XY[complete.cases(XY), ]
      y <- XY[[ yvar ]]
      X <- XY[,xvarS[[xi]]]
      full.model <- lm(y ~., data = X)
      step.model <- stepAIC(full.model, direction = "both", trace = FALSE)
      print(summary(step.model))
      
    }
  }
}
