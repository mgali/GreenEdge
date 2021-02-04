# RDA analysis of DMS and DMSPt vs. pigments and physical variables following Ch. 6 of book Numerical Ecology with R

library(RColorBrewer)
library(tidyr)
library(dplyr)
library(ade4)
library(vegan)
# install.packages("ggvegan") # package ‘ggvegan’ is not available (for R version 3.5.1). Needed for nicer RDA plotting

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
prof <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
surf <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)

# Output settings
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/RDA/"
exportimg <- T
cp_norm <- F # Not relevant additional info.

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

# Rename variables (columns)
prof <- rename(prof,
               PAR24h = par_d_p24h_ein_m_2_day_1,
               DMSPt = dmspt,
               DMS = dms,
               Temp = temp,
               Sal = sal,
               ANP = anp,
               TChl_a = tchla,
               Chl_b = chlb,
               Chl_c2 = chlc2group,
               Chl_c3 = chlc3,
               Fuco = fuco,
               Hex_fuco = hex,
               But_fuco = but,
               DD = dd,
               TCar = tcar,
               Pheo_a = phdaSUM,
               Peri = peri,
               But_fuco_like = but19_like,
               Allo = allo,
               Pras = pras.1,
               VAZ = vaz
)

# ------------------------------------------------------------------------
# Prepare matrices
# Pigments with less than 50% available data not considered

# Define subsets
yvarS <- c("DMSPt","DMS")
xvarS <- list(physics = c("Temp","Sal","N2","ANP","PAR24h"), # ,"cpsmooth1" 
              pigments = c("TChl_a","Chl_b","Chl_c2","Chl_c3","Fuco","Hex_fuco","But_fuco","DD","TCar","Pheo_a","Peri","But_fuco_like","Allo","Pras","VAZ"))
varsall <- c( "stn","depth",yvarS, xvarS[["physics"]], xvarS[["pigments"]] )


# Set pigments = NA to 0 assuming they were below detection limit. Otherwise important samples excluded
prof[, xvarS[["pigments"]] ][is.na(prof[, xvarS[["pigments"]] ]) & prof$stn %in% c(507,512,519,612,615)] <- 0

# Divide sulfur and pigments by cp (cp_norm key)
if (cp_norm) {
  cp <- prof$cpsmooth1
  NUM <- prof[, c( yvarS, xvarS[["pigments"]] ) ]
  prof[, c( yvarS, xvarS[["pigments"]] ) ] <- NUM / as.data.frame( matrix ( cp, nrow = length(cp), ncol = dim(NUM)[2] ))
}

# Select complete cases (rows)
coca <- complete.cases( prof[,varsall] )
PRE <- prof[coca,varsall]
# rownames(PRE) <- PRE$stn # duplicated row names not allowed

# Merge with profiles to get clustering coefficient and SIC classification
PRE <- merge(x = PRE, y = surf[,c("stn","OWD")], all.x = T, all.y = F, by = "stn", suffixes = "")
toremove <- c( FALSE, !diff(PRE$stn) & !diff(PRE$depth) )
PRE <- PRE[ !toremove ,]

# Add columns that will be used afterwards
PRE$owd_class = cut(PRE$OWD, breaks = c(-35,-3.5,3.5,35), labels = c("ICE","MIZ","OW")) # OWD class
PRE$depth[PRE$depth == 0.7] <- 1                                                        # paste stn and depth for labels. Set surface to 1 for simplicity
PRE$stn_depth <- paste(PRE$stn, PRE$depth, sep = "_")

# ==== Check normality ====
# Shapiro-Wilk method?

# ==== Transform and center ====
# PRE <- log10(PRE + 1)
# PRE <- as.data.frame(scale(PRE, center = T, scale = T))

# Define X and Y matrices for RDA. Transform and center if not done in previous step
Y <- as.data.frame( scale( log10(PRE[,yvarS] + 1), center = T , scale = T) )
X <- as.data.frame( scale( log10(PRE[,xvarS[["physics"]]] + 1) , center = T, scale = T) )
W <- as.data.frame( scale( log10(PRE[,xvarS[["pigments"]]] + 1) , center = T, scale = T) )

# ------------------------------------------------------------------------
# Variation partitioning: 3 RDAs in one: Y vs. X, W vs. X, and Y vs X&W.
# The common explanatory power of X&W on Y can be estimated by difference
# Main drawback is that the model is linear whereas DMS(P) may show unimodal
# rather than linear response to some physical variables (eg ANP)
# ------------------------------------------------------------------------
vp <- varpart(Y, X, W)
# print(vp)

# ------------------------------------------------------------------------
# RDA: Choose scaling. default is 2
# ------------------------------------------------------------------------
scsc <- 2
const <- 3.427439 # scaling factor to get cases (rows) projected at right place. Obtained by printing rda.sc, don't know how to access it directly

# ------------------------------------------------------------------------
###################### RDA1: Sulfur versus physics #######################
# ------------------------------------------------------------------------
rda1 <- rda(Y ~ Temp+Sal+N2+ANP+PAR24h,
            data = X) # formula call: once exploration has been done
print(rda1)
print(coef(rda1))
rda1.sc <- scores(rda1, choices=1:2, scaling=scsc, display = "sp") # Compute scores

# Compute adjusted R2
R2_1 <- RsquareAdj(rda1)$adj.r.squared
# print(R2_1)

# Variance inflation factors (VIF), which measure the proportion by which the variance of a
# regression coefficient is inflated in the presence of other explanatory variables.
VIF1 <- vif.cca(rda1)
# print(VIF1)

# Global test of the RDA result
testglob1 <- anova.cca(rda1, step = 1000)
# print(testglob1)

# Tests of all canonical axes
testbyax1 <- anova.cca(rda1, by="axis", step=1000)
# print(testbyax1)

# ------------------------------------------------------------------------
###################### RDA1: Sulfur versus pigments ######################
# ------------------------------------------------------------------------
rda2 <- rda(Y ~ TChl_a+Chl_b+Chl_c2+Chl_c3+Fuco+Hex_fuco+But_fuco+TCar+Pheo_a+Peri+But_fuco_like+Allo+Pras+VAZ+DD, # +VAZ+DD
            data = W) # formula call: once exploration has been done
# print(rda2)
# print(coef(rda2))
rda2.sc <- scores(rda2, choices=1:2, scaling=scsc, display="sp") # Compute scores

# Compute adjusted R2
R2_2 <- RsquareAdj(rda2)$adj.r.squared
# print(R2_2R2_2)

# Variance inflation factors (VIF), which measure the proportion by which the variance of a
# regression coefficient is inflated in the presence of other explanatory variables.
VIF2 <- vif.cca(rda2)
# print(VIF2)

# Global test of the RDA result
testglob2 <- anova.cca(rda2, step = 1000)
# print(testglob2)

# Tests of all canonical axes
testbyax2 <- anova.cca(rda2, by="axis", step=1000)
# print(testbyax2)


# ------------------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------------------

# Palette
pal <- colorRampPalette(brewer.pal(9, "Spectral"))
col <- list(ICE = pal(n = 21)[21],
            MIZ = pal(n = 21)[18],
            OW = pal(n = 21)[5])


if (exportimg) {
  png(filename = ifelse(cp_norm,
                        "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/RDA/plots_last/Fig_RDA_cpnorm.png",
                        "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/RDA/plots_last/Fig_RDA.png"),
      width = 16, height = 12, units = 'cm', pointsize = 9,
      bg = 'white', res = 600, type = 'cairo')
}

# Multipanel
ma <- matrix(data = 1, nrow = 4, ncol = 16)
mb <- matrix(data = 2, nrow = 8, ncol = 7)
mc <- matrix(data = 3, nrow = 8, ncol = 9)
layout(rbind(ma, cbind(mb, mc))) # layout.show(10)
par(oma = c(1,1,0.5,0.5))

# VARPART plot: Venn diagram of explained variances
par(mar = c(3,15,1,15))
plot(vp, digits=1, bg = c("darkgray", "darkslategray"), Xnames = c("Physics","Pigments"), bty = "n") 

# Triplot RDA1
par(mar = c(6,6,1,1))
plot(rda1,
     scaling=scsc,
     display = c("sp","cn"), # c("sp","lc","cn"),
     xlim = c(-1.8,1),
     ylim = c(-1, 1),
     main = paste("DMS and DMSPt vs. Physics"),
     xlab = paste0("RDA1, ",round(100*testbyax1$Variance[1]/sum(testbyax1$Variance),digits = 1),"% variance"),
     ylab = paste0("RDA2, ",round(100*testbyax1$Variance[2]/sum(testbyax1$Variance),digits = 1),"% variance")
)
arrows(0,0,rda1.sc[,1],rda1.sc[,2], length = 0, lty = 3, lwd = 0.5, col = "red")
# Samples (cases, rows)
x <- rda1$CCA$u[,1] * const
y <- rda1$CCA$u[,2] * const
for ( cc in levels(PRE$owd_class) ) {
  iplot <- PRE$owd_class == cc
  points(x[iplot], y[iplot], pch = 20, cex = 1, col = col[[cc]])
}
itext <- PRE$DMS > 20 | PRE$DMSPt > 200
text(x[itext], y[itext], PRE$stn_depth[itext], pos = 4, offset = 0.3, cex = 0.85, col = "black") # [iplot & itext]

# Triplot RDA2
par(mar = c(6,6,1,1))
plot(rda2,
     scaling=scsc,
     display = c("sp","cn"), # c("sp","lc","cn"),
     xlim = c(-2.2,0.7),
     ylim = c(-1, 1),
     main = ifelse(cp_norm,
                   paste("DMS and DMSPt vs. Pigments, Cp-normalized"),
                   paste("DMS and DMSPt vs. Pigments")),
     xlab = paste0("RDA1, ",round(100*testbyax2$Variance[1]/sum(testbyax2$Variance),digits = 1),"% variance"),
     ylab = paste0("RDA2, ",round(100*testbyax2$Variance[2]/sum(testbyax2$Variance),digits = 1),"% variance")
)
arrows(0,0,rda2.sc[,1],rda2.sc[,2], length = 0, lty = 3, lwd = 0.5, col = "red")
# Samples (cases, rows)
x <- rda2$CCA$u[,1] * const
y <- rda2$CCA$u[,2] * const
for ( cc in levels(PRE$owd_class) ) {
  iplot <- PRE$owd_class == cc
  points(x[iplot], y[iplot], pch = 20, cex = 1, col = col[[cc]])
}
itext <- PRE$DMS > 20 | PRE$DMSPt > 200
text(x[itext], y[itext], PRE$stn_depth[itext], pos = 4, offset = 0.3, cex = 0.85, col = "black") # [iplot & itext]

if (exportimg) {dev.off()}
