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

# Add ratios
prof$dms2dmspt <- prof$dms/prof$dmspt                                # dms/dmspt ratio
prof$dmspt2tchla <- prof$dmspt/prof$tchla                            # dmspt/tchla ratio
prof$dmspt2cp <- prof$dmspt/prof$cpsmooth1                           # dmspt/cp ratio
prof$chlb2tchla <- prof$chlb/prof$tchla                              # chlb to tchla
prof$chlc2_2_tchla <- prof$chlc2group/prof$tchla                     # chlc2 to tchla (chlc2 widespread pigment)
prof$chlc3_2_tchla <- prof$chlc3/prof$tchla                          # chlc3 to tchla (Phaeocystis proxy?)
prof$fuco_2_tchla <- prof$fuco/prof$tchla                            # Fucoxanthin to tchla (fuco widespread pigment)
prof$hex_2_tchla <- prof$hex/prof$tchla                              # 19-but to tchla (proxy of what? Phaeocystis? other hapto?)
prof$but_2_tchla <- prof$but/prof$tchla                              # 19-but to tchla (proxy of what? Phaeocystis? other hapto?)
prof$peri_2_tchla <- prof$peri/prof$tchla                            # Peridinin to tchla (peri in dinos)
prof$but19like_2_tchla <- prof$but19_like/prof$tchla                 # chlc3 to tchla (Phaeocystis proxy?)
prof$tcar2tchla <- prof$tcar/prof$tchla                              # beta carotene to tchla
prof$dd2tchla <- prof$dd/prof$tchla                                  # D+D xantophyll cycle pigments to tchla
prof$vaz2tchla <- prof$vaz/prof$tchla                                # VAZ xantophyll cycle pigments to tchla
prof$phdaSUM2tchla <- prof$phdaSUM/prof$tchla                        # Phaeophorbide a to TChl (HPLC)

# ------------------------------------------------------------------------
# Prepare regression and run in loop
# Pigments with less than 50% available data

yvarS <- c("dmspt","dms","dms2dmspt")
xvarS <- list(pigments = c("chlc3","chlc2group","chldaSUM","peri","phdaSUM","but","fuco","neo","pras.1","hex",
                           "dd","allo","lut","chlb","tchla","phytnaSUM","tcar","but19_like","vaz"),
              # physics = c("temp","sal","N2","cpsmooth1","anp","par_d_p24h_ein_m_2_day_1"), # ,"dmspt"
              # physics = c("temp","sal","N2","cpsmooth1","anp","par_d_p24h_ein_m_2_day_1","dmspt"))
              # physics = c("cpsmooth1","dmspt"))
              # pigments4dmsp = c("tchla","chlb","chlc2group","chlc3","fuco","hex","tcar")) # Best subset for DMSPt. No improvement adding chlc3, 19'but-like, but or peri (although the last enters with marginally significant)
              # pigments4dmspBIS = c("tchla","chlb","chlc2group","chlc3","fuco","hex","tcar","but19_like")) # Best subset for DMSPt. No improvement adding chlc3, 19'but-like, but or peri (although the last enters with marginally significant)
              pigments4dmspBIS = c("tchla","chlb","chlc2group","chlc3","fuco","hex","tcar","phdaSUM")) # Best subset for DMSPt. No improvement adding chlc3, 19'but-like, but or peri (although the last enters with marginally significant)
              # pigments4dms = c("tchla","chlb","chlc2group","fuco","hex","tcar")) # Best subset for DMS. Replacing tchla by chlb makes a difference. Adding peri: marginally significant and improves R2 to 0.82
              # pigments4dmsBIS = c("tchla","chlb","chlc2group","fuco","hex","tcar","but19_like")) # Best subset for DMS BIS: idem adding 19'but-like. Very high R2 but N drops by 1/3
              # pigments4dms2dmspt = c("chlb2tchla","chlc2_2_tchla","chlc3_2_tchla","fuco_2_tchla","hex_2_tchla","but_2_tchla","peri_2_tchla","but19like_2_tchla","tcar2tchla","dd2tchla","vaz2tchla","phdaSUM2tchla"))
              # pigments4dms2dmsptBIS = c("fuco_2_tchla","hex_2_tchla","but19like_2_tchla","tcar2tchla","phdaSUM2tchla"))

for (mm in c("spearman","pearson")) {
  for (yvar in yvarS[3]) { #[1]
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
      # XY <- as.data.frame(scale(XY, center = T, scale = T)) # Standardize X and y
      y <- XY[[ yvar ]]
      X <- XY[,xvarS[[xi]]]
      X <- as.data.frame(scale(X, center = T, scale = T)) # Standardize X only (intercept will be mean y)
      
      full.model <- lm(y ~., data = X)
      step.model <- stepAIC(full.model, direction = "both", trace = FALSE)
      print(summary(step.model))
      
    }
  }
}
