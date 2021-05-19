# RDA analysis of DMS and DMSPt vs. pigments and physical variables following Ch. 6 of book Numerical Ecology with R

# install.packages("tidyr")
# install.packages("vegan")
# install.packages("dplyr")
# install.packages("ade4")
# install.packages("packfor") # package ‘packfor’ is not available (for R version 3.5.1)
# install.packages("ellipse")
# install.packages("FactoMineR")
library(tidyr)
library(dplyr)
library(ade4)
library(vegan)
# library(packfor)
# library(MASS)
# library(ellipse)
# library(FactoMineR)

# Load data
genpath <- '~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/RDA/'
prof <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)

# Output settings
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/RDA/"
exportimg <- F

# Rename DMS variable and remove unnecessary DMS variables
prof$dms <- prof$dms_consens_cf68
toinclude <- names(prof)[grep("dms",names(prof), invert = T)]
toinclude <- c(toinclude,"dms","dmspt")
prof <- prof[,toinclude]
prof <- prof[,names(prof)[grep("cf",names(prof), invert = T)]]

# Rename PAR
prof$PAR24h <- prof$par_d_p24h_ein_m_2_day_1

# Remove data where no DMS or DMSPt are available (needed for assignment of new values in next step)
prof <- prof[(!is.na(prof$dms) | !is.na(prof$dmspt)) & !is.na(prof$depth),]

# !!!!!! Correct DMS and DMSPt in stn 519 surface !!!!!
prof[prof$stn==519 & prof$depth==0.7,c("dms","dmspt")] <- c(11.42,79.9)

# # Count proportion of empty cells in each column (mostly used to exclude some pigments or group them into functional units as DD or VAZ xanthophylls)
# noNAcount <- sapply(prof, function(x) round( sum(!is.na(x) & (!is.na(prof$dms) | !is.na(prof$dmspt))), digits = 2) )

# Group DD and VAZ cycles
prof$dd <- rowSums(prof[,c("diadino","diato")], na.rm = T)
prof$vaz <- rowSums(prof[,c("zea","anthera","viola")], na.rm = T)

# Temperature to Kelvin
prof$temp <- prof$temp + 273.15

# Negative Brunt-Vaisala to 0
prof$N2[prof$N2 < 0] <- 0

# NA 19-hex, peri, allo and pras to 0
prof$but19_like[is.na(prof$but19_like)] <- 0
prof$peri[is.na(prof$peri)] <- 0
prof$allo[is.na(prof$allo)] <- 0
prof$pras[is.na(prof$pras)] <- 0

# ------------------------------------------------------------------------
# Prepare matrices
# Pigments with less than 50% available data not considered

# yvarS <- c("dmspt")
# yvarS <- c("dms")
yvarS <- c("dmspt","dms")

xvarS <- list(physics_and_dmspt = c("temp","sal","N2","anp","PAR24h"),
              physics = c("temp","sal","N2","anp","PAR24h"),
              pigments = c("tchla","chlb","chlc2group","chlc3","fuco","hex","but","dd","tcar","phdaSUM","peri","but19_like","allo","pras","vaz"), # ,"phdaSUM"
              # pigments = c("tchla","chlb","chlc2group","chlc3","fuco","hex","but","dd","tcar","phdaSUM","allo","pras","vaz"), # ,"phdaSUM"
              pigments_major = c("tchla","chlb","chlc2group","chlc3","fuco","hex","but","dd","tcar","phdaSUM"), # ,"phdaSUM"
              pigments_reduced = c("chlb","chlc2group","hex","but","dd","tcar"), # Excluding pigments with high collinearity: VIF>20 (log psace data)
              pigments_19butlike = c("tchla","chlb","chlc2group","chlc3","fuco","hex","but","dd","tcar","but19_like"),
              pigments_19butlike_reduced = c("chlb","fuco","hex","tcar","but19_like"),
              bio = c("tchla","chlb","chlc2group","chlc3","fuco","hex","but","dd","tcar","cpsmooth1"),
              bio_19butlike = c("tchla","chlb","chlc2group","chlc3","fuco","hex","but","dd","tcar","but19_like","cpsmooth1"))
xx <- "physics"; #yvarS <- c(xvarS$bio, yvarS) # yvarS <- c(xvarS$bio)
# xx <- "pigments"
# xx <- "pigments_major"
# xx <- "pigments_reduced"
# xx <- "pigments_19butlike"
# xx <- "pigments_19butlike_reduced"
# xx <- "bio"
# xx <- "bio_19butlike"

# Select complete cases (rows)
coca <- complete.cases(prof[,c(yvarS,xvarS[[xx]])])
PRE <- prof[coca,c(yvarS,xvarS[[xx]])]

# Create unique names for each sample
# prof = prof %>% 
#   tidyr::unite(rn, stn, cast, depth, sep = "_", remove = FALSE)
# rownames(PRE) <- prof[coca,"stn"] # Not unique!
# rownames(PRE) <- 1:dim(PRE)[1]
rownames(PRE) <- NULL
  
# ==== Test for normality ====
# for (vv in names(prof)) {
#   if (sum(diff(prof[[vv]]), na.rm = T) != 0) {
#     sw <- shapiro.test(prof[[vv]]) # Shapiro-Wilk test
#   }
#   if (sw$p.value >= 0.05) {
#     print(paste0(vv," is normal"))
#   } else {
#     print(sw$p.value)
#   }
# }

# ==== Transform and center ====
PRE <- log10(PRE + 1)
PRE <- as.data.frame(scale(PRE, center = T, scale = T))
# ==============================

# Define X and Y matrices for RDA
Y <- PRE[,yvarS]
names(Y) <- c("DMSPt","DMS")

X <- PRE[,xvarS[[xx]]]
if (xx=="pigments") {
  names(X) <- c("tchla","chlb","chlc2","chlc3","fuco","hex","but","dd","tcar","pheo","peri","but_like","allo","pras","vaz")
  # names(X) <- c("tchla","chlb","chlc2","chlc3","fuco","hex","but","dd","tcar","pheo","allo","pras","vaz")
}

# Triplot scaling. Default is 2
scsc <- 2 # Choose 1 or 2

# ------------------------------------------------------------------------
# Simple, exploratory RDAs (Redundancy analysis)
# ------------------------------------------------------------------------

# # Simple RDA (no qualitative factors)
# rda1 <- rda(Y,X) # default (simple) call: useful for exploratory purposes
# print(rda1)
# print(coef(rda1))
# rd1.sc <- scores(rda1, choices=1:2, scaling=scsc, display="sp") # Compute scores
# 
# # Triplot
# plot(rda1,
#      scaling=scsc,
#      display = c("sp","lc","cn"),
#      main = paste("Triplot ",paste(yvarS,collapse = "_"),sep="_"))
# arrows(0,0,rd1.sc[,1],rd1.sc[,2], length = 0, lty = 1, lwd = 2, col = "red")
# 
# # Additional statistics
# 
# # Compute adjusted R2
# R2 <- RsquareAdj(rda1)$adj.r.squared
# print(R2)
# 
# # Variance inflation factors (VIF), which measure the proportion by which the variance of a
# # regression coefficient is inflated in the presence of other explanatory variables.
# VIF <- vif.cca(rda1)
# print(VIF)
# 
# # Global test of the RDA result
# testglob <- anova.cca(rda1)

# ------------------------------------------------------------------------
# More targetted RDA: Sulfur versus pigments
# ------------------------------------------------------------------------
# rda2 <- rda(Y ~ tchla+chlb+chlc2group+chlc3+fuco+hex+but+dd+tcar+phdaSUM,
#             data = X) # formula call: once exploration has been done
# print(rda2)
# print(coef(rda2))
# rda2.sc <- scores(rda2, choices=1:2, scaling=scsc, display="sp") # Compute scores
# 
# # Compute adjusted R2
# R2 <- RsquareAdj(rda2)$adj.r.squared
# print(R2)
# 
# # Variance inflation factors (VIF), which measure the proportion by which the variance of a
# # regression coefficient is inflated in the presence of other explanatory variables.
# VIF <- vif.cca(rda2)
# print(VIF)
# 
# # Global test of the RDA result
# testglob <- anova.cca(rda2, step = 1000)
# print(testglob)
# 
# # Tests of all canonical axes
# testbyax <- anova.cca(rda2, by="axis", step=1000)
# print(testbyax)
# 
# # Triplot
# plot(rda2,
#      scaling=scsc,
#      display = c("sp","lc","cn"),
#      main = paste("Triplot ",paste(yvarS,collapse = "_"),sep="_"))
# arrows(0,0,rda2.sc[,1],rda2.sc[,2], length = 0, lty = 3, lwd = 2, col = "red")


# ------------------------------------------------------------------------
# More targetted RDA
# ------------------------------------------------------------------------
rda3 <- rda(Y ~ temp+sal+N2+anp+PAR24h,
# rda3 <- rda(Y ~ tchla+chlb+chlc2+chlc3+fuco+hex+but+dd+tcar+pheo+allo+pras+vaz,
# rda3 <- rda(Y ~ tchla+chlb+chlc2+chlc3+fuco+hex+but+dd+tcar+pheo+allo+pras+but_like+peri+vaz,
# rda3 <- rda(Y ~ tchla+chlb+chlc2group+chlc3+fuco+hex+but+dd+tcar+phdaSUM,
            data = X) # formula call: once exploration has been done
print(rda3)
print(coef(rda3))
rda3.sc <- scores(rda3, choices=1:2, scaling=scsc, display="sp") # Compute scores

# Compute adjusted R2
R2 <- RsquareAdj(rda3)$adj.r.squared
print(R2)

# Variance inflation factors (VIF), which measure the proportion by which the variance of a
# regression coefficient is inflated in the presence of other explanatory variables.
VIF <- vif.cca(rda3)
print(VIF)

# Global test of the RDA result
testglob <- anova.cca(rda3, step = 1000)
print(testglob)

# Tests of all canonical axes
testbyax <- anova.cca(rda3, by="axis", step=1000)
print(testbyax)

# Triplot
scsc <- 2 # Choose scaling. default is 2
plot(rda3,
     scaling=scsc,
     pch = 2,
     display = c("sp","lc","cn"),
     # main = paste("Triplot ",paste(yvarS,collapse = "_"),sep="_"),
     # main = paste("DMS and DMSPt vs. pigments"),
     main = paste("DMS and DMSPt vs. physics"),
     xlab = paste0(round(100*testbyax$Variance[1]/sum(testbyax$Variance),digits = 1),"% variance"),
     ylab = paste0(round(100*testbyax$Variance[2]/sum(testbyax$Variance),digits = 1),"% variance"),
     )
arrows(0,0,rda3.sc[,1],rda3.sc[,2], length = 0, lty = 3, lwd = 1, col = "red")
