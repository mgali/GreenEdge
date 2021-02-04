# MAKE PLOTS OF DMS AND DMSP VS. PHYTO COUNTS AND ADDITIONAL Z VARIABLES

library(RColorBrewer)
library(dplyr)

# Load data
genpath <- '~/Desktop/GreenEdge/GCMS/'
prof.all <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), header = T)
surf.all <- read.csv(file = paste0(genpath,'GE2016.casts.ALLSURF.csv'), header = T)

# Exporting image?
exportimg <- T
xrange <- c(0.6,-0.6)
opath <- "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_ salinity_corr/"

# ---------------------
pal <- colorRampPalette(brewer.pal(12, "Paired"))              # Color palette (Consistent with Matlab plots)
col <- pal(n = 12)
plotres <- 600                                                  # resolution dpi

# ---------------------
# Rename DMS variable and remove unnecessary DMS variables
prof.all$dms <- prof.all$dms_consens_cf68
toinclude <- names(prof.all)[grep("dms",names(prof.all), invert = T)]
toinclude <- c(toinclude,"dms","dmspt")
prof.all <- prof.all[,toinclude]

# Remove data where no DMS or DMSPt are available
prof.all <- prof.all[(!is.na(prof.all$dms) | !is.na(prof.all$dmspt)) & !is.na(prof.all$depth),]

# !!!!!! Correct DMS and DMSPt in stn 519 surface !!!!!
# prof.all[prof.all$stn==519 & prof.all$depth==0.7,c("dms","dmspt")] <- c(3.93,79.9)
prof.all[prof.all$stn==519 & prof.all$depth==0.7,c("dms","dmspt")] <- c(11.42,79.9)

# Add MIZ classification
icecrit1 <- 0.15
icecrit2 <- 0.70
ICE <- surf.all[,c("SICm2d","SICm1d","SICday")]
icemin <- apply(ICE, 1, min, na.rm = T) # Min-max SIC criterion
icemax <- apply(ICE, 1, max, na.rm = T) # Min-max SIC criterion
icemean <- apply(ICE, 1, min, na.rm = T) # Mean concentration criterion
surf.all$sic_class <- NA
surf.all$sic_class[icemax<icecrit1] <- "OW"
surf.all$sic_class[icemin>icecrit2] <- "ICE"
surf.all$sic_class[icemax>=icecrit1 & icemin<=icecrit2] <- "MIZ"

# Merge with profiles to get clustering coefficient and SIC classification
pplot <- merge(x = prof.all, y = surf.all, all.x = T, all.y = F, by = 'stn', suffixes = "")

# Remove duplicated columns with NA in their names
pplot <- pplot[,grep("NA",names(pplot), invert = T)]

# Remove duplicated rows 
dd <- duplicated(pplot[,c("dmspt","dms","cast","depth")]) | is.na(pplot$dmspt)
pplot <- pplot[!dd,]

# # Hide data from stn 400 (either entire or just surface)
# # pplot[pplot$stn<=400,] <- NA
# pplot[pplot$stn<=400 & pplot$depth < 5,] <- NA

# Change units of N2 from s-2 to h-1
pplot$N2 <- sqrt(pplot$N2) * 3600

# Fill NA with zeros in pigments frequently below LOD and only if CHl is not NA
pplot[!is.na(pplot$tchla),c("chlc3","but19_like","peri")][is.na(pplot[!is.na(pplot$tchla),c("chlc3","but19_like","peri")])] <- 0

# ---------------------
# Bin profiles by station categories
df2bin <- pplot
z_class <- cut(df2bin$depth, breaks = c(0,9,21,41,81), labels = c(0,1,2,3))
st_class <- list(sic_class = pplot$sic_class,
                 owd_class = cut(df2bin$OWD, breaks = c(-35,-3.5,3.5,35), labels = c("ICE","MIZ","OW")))

# # Uncomment to select MIZ data only
# pplot <- pplot[st_class$sic_class %in% c("MIZ"),]
# # pplot <- pplot[pplot$OWD>-5 & pplot$OWD<5,]
# xrange <- c(0.8,-0.8)

# ---------------------
# Plot settings. EXPLORATORY
yvarS <- list(dmspt = "DMSPt",
              dms = "DMS",
              chlc2group = "Chl c2",
              tchla = "TChl a",
              but19_like = "But-fuco-like",
              chlc3 = "Chl c3",
              but = "But-fuco",
              hex = "Hex-fuco",
              fuco = "Fuco",
              peri = "Peri",
              phdaSUM = "Pheo a",
              cpsmooth1 = "Cp"
)

# ---------------------
# Correlation between salinity and yvarS
maxdepth <- c(41, 21, 10, 5)
names(col) <- names(yvarS)
selvars <- c("dms","dmspt","tchla","cpsmooth1","chlc3","but19_like")
yvarS <- yvarS[names(yvarS) %in% selvars]
col <- col[names(col) %in% selvars]

CC <- list()

if (exportimg) {png(filename = paste0(opath,"Fig_salinity_corr.png"), width = 8, height = 6, units = 'cm', pointsize = 6, bg = 'white', res = plotres, type = 'cairo')}

par(mar=c(4,4,3,8))

pplot <- pplot[complete.cases(pplot[,selvars]),]

for (y in names(yvarS)) {
  
  print(y)
  C <- sapply(maxdepth, function(x) {
    iz <- pplot$depth <= x
    return(cor.test(pplot$sal[iz], pplot[iz,y], method = "spearman"))
    # return(cor.test(pplot$sal[iz], pplot[iz,y], method = "pearson"))
  })
  CC[[y]] <- cbind(C["estimate",],C["p.value",])
  
  if (y==names(yvarS)[1]) {
    plot(C["estimate",],
         maxdepth,
         cex = 1.3,
         ylim = rev(range(maxdepth)),
         xlim = xrange,
         ylab = "Maximum depth (m)",
         xlab = expression("r"[S]*" with salinity"),
         pch = 1, col=col[[y]])
    abline(v = 0, col = "gray", lty = 2)
  }
  
  isig <- C["p.value",] < 0.05
  points(C["estimate",], maxdepth, pch = 1, cex = 1.3, col=col[[y]])
  lines(C["estimate",], maxdepth, col=col[[y]])
  points(C["estimate",isig], maxdepth[isig], pch = 16, cex = 1.3, col=col[[y]])
  
}


mtext(text = "Enhanced by low salinity", side = 3, line = 1.5, at = min(xrange)/2 - 0.05)
mtext(text = "Enhanced by high salinity", side = 3, line = 1.5, at = max(xrange)/2 + 0.05)

par(xpd=TRUE)
arrows(-0.1,2,min(xrange),2, length = 0.03, angle = 30, col = "black")
arrows(+0.1,2,max(xrange),2, length = 0.03, angle = 30, col = "black")

legend(x = min(xrange)-0.05,
       y = 15,
       legend = yvarS, bty = "n",
       pch = 16,
       cex = 1.1,
       col = col)

if (exportimg) {dev.off()}
