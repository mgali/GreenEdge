# ADD PHAEOPHORBIDE TO PIGMENTS IN PROFILES DATABASE

# install.packages("readxl")
library("readxl")

# Load profiles database
genpath <- '~/Desktop/GreenEdge/GCMS/'
prof <- read.csv(file = paste0(genpath,'GE2016.profiles.ALL.OK.backup2.csv'), header = T)

# Load phaeophorbide data
pigpath <- '~/Desktop/GreenEdge/'
pigments <- readxl::read_excel(path = paste0(pigpath,"GreenEdge-Amundsen-pigments-180131.xlsx"), sheet = "Discrete", col_names = T, na = "LOD")
ppho <- data.frame(stn = pigments$Station,
                   depth = pigments$`Depth (m)`,
                   phbda = pigments$`sum Phbd a`)

# Modify station names to match numeric ones in main database
ppho$stn <- as.numeric(gsub(pattern = "G", replacement = "", ppho$stn))

# Add in loop to ensure match within plus/minus 1 m
st <- unique(ppho$stn)[unique(ppho$stn) %in% unique(prof$stn)]
prof$phbda <- 0 * prof$tchla

for (s in st) {
  
  tmp <- ppho[ppho$stn==s & !is.na(ppho$stn),]
  
  for (z in unique(prof$depth[prof$stn==s])) {
    iz <- which(prof$depth==z & prof$stn==s & !is.na(prof$depth)  & !is.na(prof$stn))
    tmpmean <- mean(tmp$phbda[(tmp$depth >= z-1) & (tmp$depth <= z+1)], na.rm = T)
    prof$phbda[iz] <- tmpmean
  }
}

# NaN to NA
prof[is.nan(prof)] <- NA

# Write data out
write.csv(x = prof, file = paste0(genpath,'GE2016.profiles.ALL.OK.csv'), row.names = F)

# ---------------
# # Checks
# plot(prof$phbda, prof$depth)
# plot(prof$tchla,prof$phbda)
# plot(prof$phbda/prof$tchla, prof$depth)



