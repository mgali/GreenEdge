# ADD ALL PIGMENTS (AGAIN) TO PROFILES DATABASE. SOME MISSING VALUES DETECTED IN OLD MATCH

# install.packages("readxl")
library("readxl")
library(dplyr)

# Load profiles database
genpath <- '~/Desktop/GreenEdge/GCMS/'
prof <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.backup3.csv'), header = T)

# Load phaeophorbide data
pigpath <- '~/Desktop/GreenEdge/'
pigments <- readxl::read_excel(path = paste0(pigpath,"GreenEdge-Amundsen-pigments-180131.xlsx"), sheet = "Discrete", col_names = T, na = "LOD")
addpig <- data.frame(cast = pigments$CTD,
                     depth = pigments$`Depth (m)`,
                     chlc3 = pigments$`Chlorophyll c3`,
                     chlc2group = pigments$`Chlorophyll c1+c2+MgDVP`,
                     chldaSUM = pigments$`sum Chld a`,
                     peri = pigments$Peridinin,
                     phdaSUM = pigments$`sum Phbd a`,
                     uriol = pigments$Uriolide,
                     but = pigments$`19'-Butanoyloxyfucoxanthin`,
                     fuco = pigments$Fucoxanthin,
                     neo = pigments$Neoxanthin,
                     pras.1 = pigments$Prasinoxanthin,
                     viola = pigments$Violaxanthin,
                     hex = pigments$`19'-Hexanoyloxyfucoxanthin`,
                     asta = pigments$Astaxanthin,
                     micral = pigments$Micromonal,
                     diadino = pigments$Diadinoxanthin,
                     anthera = pigments$Antheraxanthin,
                     allo = pigments$Alloxanthin,
                     diato = pigments$Diatoxanthin,
                     zea = pigments$Zeaxanthin,
                     lut = pigments$Lutein,
                     bchla = pigments$`Bacteriochlorophyll a`,
                     chlb = pigments$`Chlorophyll b`,
                     dvchla = pigments$`Divinyl Chlorophyll a`,
                     chla = pigments$`Chlorophyll a`,
                     tchla = pigments$`Total Chlorophyll a`,
                     phytnaSUM = pigments$`Sum Phaeophytin a`,
                     caro_like_Prasi = pigments$`Carotene-like (prasinophyte)`,
                     tcar = pigments$`Sum a+b Carotenes`,
                     but19_like = pigments$`Sum 19BF-like`,
                     hex19_likeSUM = pigments$`Sum 19HF-like`)

# # Modify station names to match numeric ones in main database
# addpig$stn <- as.numeric(gsub(pattern = "G", replacement = "", addpig$stn))

# Modify cast names to match numeric ones in main database
addpig$cast <- as.numeric(gsub(pattern = "1601", replacement = "", addpig$cast))

# Add in loop to ensure match within plus/minus 1 m
st <- unique(addpig$cast)[unique(addpig$cast) %in% unique(prof$cast)]

# Define matching columns (excluding stn and depth)
icols <- which(names(prof) %in% names(addpig))
icols <- icols[3:length(icols)]

for (s in st) {
  
  tmp <- addpig[addpig$cast==s & !is.na(addpig$cast),]
  
  zz <- unique(prof$depth[prof$cast==s & !is.na(prof$depth)])
  
  for (z in zz) {
    
    iz <- which(prof$depth==z & prof$cast==s & !is.na(prof$depth)  & !is.na(prof$cast))
    tmpmean <- tmp[(tmp$depth >= z-0.6) & (tmp$depth <= z+0.6),3:dim(tmp)[2]]
    tmpmean <- colMeans(tmpmean, na.rm = T)
    
    for (iiz in iz) {
      prof[iiz,icols] <- tmpmean
    } 
  }
}

# NaN to NA
prof <- as.matrix(prof)
prof[is.nan(prof)] <- NA
prof <- as.data.frame(prof)

# Write data out
write.csv(x = prof, file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), row.names = F)

# ---------------
# # Checks
# plot(prof$phbda, prof$depth)
# plot(prof$tchla,prof$phbda)
# plot(prof$phbda/prof$tchla, prof$depth)



